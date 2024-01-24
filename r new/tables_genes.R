library(DESeq2)
library(EnsDb.Hsapiens.v86)
library(dplyr)
library(biomaRt)
library(ggplot2)
library(readxl)
library(tidyverse)

dir.path <- "/Volumes/AG_Landthaler/amanukyan/Projects/RiboSeq/Pipeline/run/sofia/output/FeatureCounts"
getwd()

# human data 
list_files <- list.files(dir.path)
list_files <- list_files[!grepl("summary", list_files)]
rnadata <- read.table(file="/Volumes/AG_Landthaler/amanukyan/Projects/RiboSeq/Pipeline/run/sofia_rna/output/count.txt", header = T, sep = "\t")
colnames(rnadata) <- lapply(colnames(rnadata), function(x) strsplit(x, split = "_S[0-9]+_")[[1]][1])
rownames(rnadata) <- rnadata[,1]
rnadata <- rnadata[,-1]
rnadata <- rnadata[,c(6:8,14:16)]
colnames(rnadata) <- gsub("output.bam_files.", "", colnames(rnadata))

# remove low count genes
count_gene <- apply(rnadata, 1, max) 
rnadata <- rnadata[which(count_gene > 10),]
rnaCond <- c("cyto","cyto","cyto","ER","ER","ER")
rnaCond<-factor(rnaCond)
# get regular deseq2 results
dds <- DESeqDataSetFromMatrix(countData = rnadata,
                              DataFrame(rnaCond),
                              design= ~ rnaCond)
dds <- DESeq(dds)
rna_results <- as.data.frame(results(dds))
rna_results$padj[is.na(rna_results$padj)] <- 1
# get gene lengths
exons <- exonsBy(EnsDb.Hsapiens.v86, by="gene")
exons <- GenomicRanges::reduce(exons)
len <- sum(width(exons))
len <- data.frame(genes = names(len), len = len)
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
G_list <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id", "hgnc_symbol"), values=len$genes, mart= mart)
len <- len %>% 
  left_join(G_list, by = c("genes" = "ensembl_gene_id")) %>% 
  na.omit() %>% 
  dplyr::filter(hgnc_symbol != "") %>% 
  group_by(hgnc_symbol) %>%
  summarize(len = mean(len))
# rna tpm
len_rna <- len[match(rownames(rnadata), len$hgnc_symbol),]
len_rna <- na.omit(len_rna)
rnadata_tpm <- rnadata[intersect(rownames(rnadata), len_rna$hgnc_symbol),]
rnadata_tpm <- sweep(rnadata_tpm * 1000, 1, len_rna$len, FUN = "/")
rnadata_tpm <- sweep(rnadata_tpm * 1000000, 2, colSums(rnadata_tpm), FUN = "/")
# calculate log RNA mean for each condition
cond_names <- rnaCond
Mean_RNA_Data <- aggregate(t(rnadata_tpm), list(cond_names), mean)
Mean_RNA_Data <- t(Mean_RNA_Data)
colnames(Mean_RNA_Data) <- Mean_RNA_Data[1,]
Mean_RNA_Data <- data.frame(apply(Mean_RNA_Data[-1,],2,as.numeric), row.names = rownames(Mean_RNA_Data[-1,]))
LogMean_RNA_Data <- log(Mean_RNA_Data + 0.001)
rna_results$gene <- rownames(rna_results)
datax_CE <- data.frame(log2FoldChange = rna_results$log2FoldChange, Sig_RNA = as.numeric(rna_results$padj), gene = rna_results$gene)
datax_CE$Sig_RNA[is.na(datax_CE$Sig_RNA)] <- 1

#datax[datax$gene %in% C_genes & datax$log2FoldChange > -0.5, "localisation"] <- "cytosolic"
#datax[datax$gene %in% E_genes  & datax$log2FoldChange > 1, "localisation"] <- "membrane"
#datax[!(datax$gene %in% C_genes) & !(datax$gene %in% M_genes), "localisation"] <- "undefined"

datax_CE[datax_CE$log2FoldChange < 0.5, "localisation"] <- "cytosolic"
datax_CE[datax_CE$log2FoldChange > 1.5, "localisation"] <- "membrane"
datax_CE[is.na(datax_CE$localisation) | 
           !(datax_CE$localisation == "cytosolic" | datax_CE$localisation == "membrane"), "localisation"] <- "undefined"

C_E <- read_excel("/Users/sgarcia/Downloads/cyto_ER_genes.xls",sheet = 1,col_names = T, col_types = "text")
C_genes <- C_E$Symbol[C_E$localization=="cytosolic"]
E_genes <- C_E$Symbol[C_E$localization=="membrane"]

#datax_cyto_genes <- rna_results[rna_results$log2FoldChange < -0.5,]
#datax <- rna_results[rna_results$log2FoldChange > 1,]

######## HISTOGRAM PLOT ER CYTO ############


ggplot(datax_CE, aes(log2FoldChange, fill=localisation))+geom_histogram(bins=250)+
  scale_fill_manual(values=c("dodgerblue4","orange2","grey"))

cytosolic_genes <- datax_CE[datax_CE$localisation == "cytosolic",]
membrane_genes <- datax_CE[datax_CE$localisation=="membrane",]

###### addition of byotype ###########

membrane_genes_list <- membrane_genes$gene 
gene_info_membrane <- getBM(
  attributes = c("external_gene_name", "gene_biotype"),
  filters = "external_gene_name",
  values = membrane_genes_list,
  mart = ensembl
)
membrane_genes <- membrane_genes %>%
  left_join(gene_info_membrane, by = c("gene" = "external_gene_name"))%>%
  mutate(gene_biotype = coalesce(gene_biotype, "noinfo"))
membrane_genes<-membrane_genes[membrane_genes$gene_biotype=="protein_coding",]
membrane_genes <- membrane_genes %>%
  select(gene)
write.table(membrane_genes, file = "membrane_genes.txt", sep = "\t", row.names = FALSE,col.names=F,quote=F)

cytosolic_genes_list <- cytosolic_genes$gene 
gene_info_cytosol <- getBM(
  attributes = c("external_gene_name", "gene_biotype"),
  filters = "external_gene_name",
  values = cytosolic_genes_list,
  mart = ensembl
)
cytosolic_genes <- cytosolic_genes %>%
  left_join(gene_info_cytosol, by = c("gene" = "external_gene_name"))%>%
  mutate(gene_biotype = coalesce(gene_biotype, "noinfo"))
cytosolic_genes<-cytosolic_genes[cytosolic_genes$gene_biotype=="protein_coding",]
cytosolic_genes <- cytosolic_genes %>%
  select(gene)
write.table(cytosolic_genes, file = "cytosolic_genes.txt", sep = "\t", row.names = FALSE,col.names=F,quote=F)

list_genesSP<-read.table(file="/Volumes/AG_Landthaler/sofia/riboSRP72/sp_list.txt", header = T, sep = "\t")
SP_genes<-list_genesSP$Symbol
write.table(SP_genes, file = "SP_genes.txt", sep = "\t", row.names = FALSE,col.names=F,quote=F)
SRP_dependent<-c(list_genesTM$Symbol,list_genesSP)

###### comparison genes SRP-dependent and independent proteins ##########

SRP_dependent <- read_excel("/Users/sgarcia/Documents/genes_clasiffication tables/SRP-dependent_proteins.xls",sheet = 1,col_names = T, col_types = "text",skip=3)
SRP_independent <- read_excel("/Users/sgarcia/Documents/genes_clasiffication tables/SRP-independent_proteins.xls",sheet = 1,col_names = T, col_types = "text",skip=4)

SRP_depedent_common_membrane<-SRP_dependent[SRP_dependent$`Gene name` %in% membrane_genes$gene,]
SRP_depedent_common_cytosol<-SRP_dependent[SRP_dependent$`Gene name` %in% cytosolic_genes$gene,]

SRP_indepedent_common_membrane<-SRP_independent[SRP_independent$`Gene Name` %in% membrane_genes$gene,]
SRP_indepedent_common_cytosol<-SRP_independent[SRP_independent$`Gene Name` %in% cytosolic_genes$gene,]

SRP_dependent <- SRP_dependent$`Gene name`
SRP_independent <- SRP_independent$`Gene Name`
  
# Connect to BioMart and fetch the orthologs
ensembl <- useMart("ensembl")
yeast <- useDataset("scerevisiae_gene_ensembl", mart = ensembl)
human <- useDataset("hsapiens_gene_ensembl", mart = ensembl)
orthologs <- biomaRt::getLDS(attributes = c("ensembl_gene_id", "external_gene_name"), 
                      filters = "external_gene_name", 
                      values = SRP_dependent, 
                      mart = yeast, 
                      attributesL = c("hsapiens_paralog_ensembl_gene", "hsapiens_paralog_associated_gene_name"), 
                      martL = human)
  # Save to a text file
write.table(orthologs, file="human_orthologs.txt", sep="\t", row.names=FALSE, quote=FALSE)



### vector of genes down regulated in riboseq and rnaseq under srp depletion #####
data_down <- c(datax_8_down$gene, datax_ribo_8_down$gene)
######## 
library(seqinr)

# Set up the Ensembl mart
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# List of gene names (example)
genes <- c("BRCA1", "TP53", "EGFR")

# Fetch protein sequences
sequences <- getBM(attributes = c("hgnc_symbol", "peptide"), 
                   filters = "hgnc_symbol", 
                   values = genes, 
                   mart = ensembl)
# Kyte-Doolittle scale
kd_scale <- c(A=1.8, R=-4.5, N=-3.5, D=-3.5, C=2.5, Q=-3.5, E=-3.5, G=-0.4, H=-3.2, I=4.5, L=3.8, K=-3.9, M=1.9, F=2.8, P=-1.6, S=-0.8, T=-0.7, W=-0.9, Y=-1.3, V=4.2)

# Sliding window analysis using the Kyte-Doolittle scale
for (i in 1:nrow(sequences)) {
  cat(paste("Hydropathy plot for:", sequences$hgnc_symbol[i]), "\n")
  hydropathy <- window(as.character(sequences$peptide[i]), w = window_size, sc = kd_scale, plot.it = TRUE, xlab = "Residue Number", ylab = "Hydropathy", main = paste("Kyte-Doolittle Hydropathy Plot for", sequences$hgnc_symbol[i]))
  dev.flush()
  Sys.sleep(5)  # Pause for 5 seconds
}


# Manual calculation of Kyte-Doolittle hydropathy
computeHydropathy <- function(sequence, window_size, kd_scale) {
  n <- nchar(sequence)
  scores <- numeric(n - window_size + 1)
  for (i in 1:(n - window_size + 1)) {
    subseq <- substr(sequence, i, i + window_size - 1)
    scores[i] <- mean(sapply(strsplit(subseq, NULL)[[1]], function(aa) kd_scale[aa]))
  }
  return(scores)
}

# Compute hydropathy for the sequence
hydropathy_scores <- computeHydropathy(test_seq, window_size, kd_scale)

# Plot the scores
plot(hydropathy_scores, type="l", xlab="Residue Number", ylab="Hydropathy", main="Kyte-Doolittle Hydropathy Plot for BRCA1", ylim=c(min(kd_scale), max(kd_scale)))