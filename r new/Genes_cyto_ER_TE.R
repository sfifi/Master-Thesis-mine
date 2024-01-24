library(DESeq2)
library(EnsDb.Hsapiens.v86)
library(dplyr)
library(biomaRt)
library(ggplot2)
library(readxl)


# library(riborex)
source("/Volumes/AG_Landthaler/sofia/DESeq2Rex.R")


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


library(readxl)
library(tidyverse)

C_E <- read_excel("/Users/sgarcia/Downloads/cyto_ER_genes.xls",sheet = 1,col_names = T, col_types = "text")
C_genes <- C_E$Symbol[C_E$localization=="cytosolic"]
E_genes <- C_E$Symbol[C_E$localization=="membrane"]

#datax[datax$gene %in% C_genes & datax$log2FoldChange > -0.5, "localisation"] <- "cytosolic"
#datax[datax$gene %in% E_genes  & datax$log2FoldChange > 1, "localisation"] <- "membrane"
#datax[!(datax$gene %in% C_genes) & !(datax$gene %in% M_genes), "localisation"] <- "undefined"

datax_CE[datax_CE$log2FoldChange < 0.5, "localisation"] <- "cytosolic"
datax_CE[datax_CE$log2FoldChange > 1.5, "localisation"] <- "membrane"
datax_CE[is.na(datax_CE$localisation) | 
           !(datax_CE$localisation == "cytosolic" | datax_CE$localisation == "membrane"), "localisation"] <- "undefined"


#datax_cyto_genes <- rna_results[rna_results$log2FoldChange < -0.5,]
#datax <- rna_results[rna_results$log2FoldChange > 1,]

######## HISTOGRAM PLOT ER CYTO ############


ggplot(datax_CE, aes(log2FoldChange, fill=localisation))+geom_histogram(bins=250)+
  scale_fill_manual(values=c("dodgerblue4","orange2","grey"))

cytosolic_genes <- datax_CE[datax_CE$localisation == "cytosolic",]
cytosolic_genes <- cytosolic_genes$gene
membrane_genes <- datax_CE[datax_CE$localisation=="membrane",]
membrane_genes <- membrane_genes$gene
undefined_genes <- datax_CE[datax_CE$localisation=="undefined",]
undefined_genes <- undefined_genes$gene

# human data 
list_files <- list.files(dir.path)
list_files <- list_files[!grepl("summary", list_files)]

ribodata <- read.table(file="/Volumes/AG_Landthaler/amanukyan/Projects/RiboSeq/Pipeline/run/sofia/output/FeatureCounts/cds_counts.txt", header = T, sep = "\t")
colnames(ribodata) <- lapply(colnames(ribodata), function(x) strsplit(x, split = "_S[0-9]+_")[[1]][1])
rownames(ribodata) <- ribodata[,1]
ribodata <- ribodata[,-1]
ribodata <- ribodata[,c(13:14,19,15,17,18)]
colnames(ribodata) <- gsub("output.bam_files.", "", colnames(ribodata))
colnames(ribodata) <- gsub(".bam", "", colnames(ribodata))

rnadata <- read.table(file="/Volumes/AG_Landthaler/amanukyan/Projects/RiboSeq/Pipeline/run/sofia_rna/output/count.txt", header = T, sep = "\t")
colnames(rnadata) <- lapply(colnames(rnadata), function(x) strsplit(x, split = "_S[0-9]+_")[[1]][1])
rownames(rnadata) <- rnadata[,1]
rnadata <- rnadata[,-1]
rnadata <- rnadata[,c(14,16,17:18,20:21)]
colnames(rnadata) <- gsub("output.bam_files.", "", colnames(rnadata))

####
# DE analysis ####
####

# remove low count genes
count_gene <- apply(rnadata, 1, max) 
rnadata <- rnadata[which(count_gene > 10),]
count_gene <- apply(ribodata, 1, max) 
ribodata <- ribodata[which(count_gene > 10),]
#Change names to intersect
rnadata <- rnadata[intersect(rownames(ribodata), rownames(rnadata)),]
ribodata <- ribodata[intersect(rownames(ribodata), rownames(rnadata)),]
# control treatment
riboCond <- c("ctrl","ctrl","time 4 hours","time 4 hours","time 8 hours","time 8 hours")
rnaCond <- c("ctrl","ctrl","time 4 hours","time 4 hours","time 8 hours","time 8 hours")
rnaCond<-factor(rnaCond)
riboCond<-factor(riboCond)
# get riborex data
print(head(rnadata))
print(head(ribodata))

# res.deseq2 <- riborex(rnadata, ribodata, rnaCond, riboCond, contrast = c("EXTRA3", "time 4 hours", "ctrl"))
#calculation of foldchange in TE among conditions.
dds.deseq2 <- DESeq2Rex(rnadata, ribodata, rnaCond, riboCond, contrast = c("EXTRA1", "time 4 hours", "ctrl"))
unique_riboCond <- as.character(rev(unique(riboCond)))
ribo_results <- NULL
for(i in 1:(length(unique_riboCond)-1)){
  for(j in (i+1):(length(unique_riboCond))){
    comparison <- c("EXTRA1", unique_riboCond[i], unique_riboCond[j])
    print(comparison)
    cur_results <- as.data.frame(results(dds.deseq2$dds, contrast = comparison))
    ribo_results <- rbind(ribo_results, data.frame(cur_results, gene = rownames(cur_results), comparison = paste0(comparison[2:3], collapse = "_vs_")))
  }
}

ribo_results$padj[is.na(ribo_results$padj)] <- 1
ribo_results_8<-filter(ribo_results, comparison == "time 8 hours_vs_ctrl")
ribo_results_4<-filter(ribo_results, comparison == "time 4 hours_vs_ctrl")


#Differential expression of RNA on ER among conditions
dds <- DESeqDataSetFromMatrix(countData = rnadata,
                              DataFrame(rnaCond),
                              design= ~ rnaCond)
dds <- DESeq(dds)
unique_rnaCond <- as.character(rev(unique(rnaCond)))
rna_results <- NULL
for(i in 1:(length(unique_riboCond)-1)){
  for(j in (i+1):(length(unique_riboCond))){
    comparison <- c("rnaCond", unique_rnaCond[i], unique_rnaCond[j])
    print(comparison)
    cur_results <- as.data.frame(results(dds, contrast = comparison))
    rna_results <- rbind(rna_results, data.frame(cur_results, gene = rownames(cur_results), comparison = paste0(comparison[2:3], collapse = "_vs_")))
  }
}
rna_results$padj[is.na(rna_results$padj)] <- 1
rna_results_8<-filter(rna_results, comparison == "time 8 hours_vs_ctrl")
rna_results_4<-filter(rna_results, comparison == "time 4 hours_vs_ctrl")


####

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

# ribo tpm
len_ribo <- len[match(rownames(ribodata), len$hgnc_symbol),]
len_ribo <- na.omit(len_ribo)
ribodata_tpm <- ribodata[intersect(rownames(ribodata), len_ribo$hgnc_symbol),]
ribodata_tpm <- sweep(ribodata_tpm * 1000, 1, len_ribo$len, FUN = "/")
ribodata_tpm <- sweep(ribodata_tpm * 1000000, 2, colSums(ribodata_tpm), FUN = "/")

# rna tpm
len_rna <- len[match(rownames(rnadata), len$hgnc_symbol),]
len_rna <- na.omit(len_rna)
rnadata_tpm <- rnadata[intersect(rownames(rnadata), len_rna$hgnc_symbol),]
rnadata_tpm <- sweep(rnadata_tpm * 1000, 1, len_rna$len, FUN = "/")
rnadata_tpm <- sweep(rnadata_tpm * 1000000, 2, colSums(rnadata_tpm), FUN = "/")

# calculate log TE mean for each condition
cond_names <- riboCond
common_genes <- intersect(rownames(ribodata_tpm), rownames(rnadata_tpm))
TE_Data <- ribodata_tpm[common_genes,]/rnadata_tpm[common_genes,]
Mean_TE_Data <- aggregate(t(TE_Data), list(cond_names), mean)
Mean_TE_Data <- t(Mean_TE_Data)
colnames(Mean_TE_Data) <- Mean_TE_Data[1,]
Mean_TE_Data <- data.frame(apply(Mean_TE_Data[-1,],2,as.numeric), row.names = rownames(Mean_TE_Data[-1,]))
LogMean_TE_Data <- log(Mean_TE_Data + 0.001)

# calculate log RNA mean for each condition
cond_names <- rnaCond
Mean_RNA_Data <- aggregate(t(rnadata_tpm), list(cond_names), mean)
Mean_RNA_Data <- t(Mean_RNA_Data)
colnames(Mean_RNA_Data) <- Mean_RNA_Data[1,]
Mean_RNA_Data <- data.frame(apply(Mean_RNA_Data[-1,],2,as.numeric), row.names = rownames(Mean_RNA_Data[-1,]))
LogMean_RNA_Data <- log(Mean_RNA_Data + 0.001)

datax_ribo_8 <- data.frame(LogFC_RNA = ribo_results_8$log2FoldChange, Sig_RNA = as.numeric(ribo_results_8$padj), gene = ribo_results_8$gene)
datax_ribo_8$Sig_RNA[is.na(datax_ribo_8$Sig_RNA)] <- 1
temp_8 <- rep("NonSig", nrow(datax_ribo_8))
temp_8[datax_ribo_8$Sig_RNA < 0.05 & abs(datax_ribo_8$LogFC_RNA) > 1] <- "Sig"
datax_ribo_8$Sig <- factor(temp_8, levels = c("NonSig", "Sig"))

datax_ribo_8$diffexpressed <- "NO"
datax_ribo_8$diffexpressed[datax_ribo_8$LogFC_RNA > 1 & datax_ribo_8$Sig == "Sig"] <- "UP"
datax_ribo_8$diffexpressed[datax_ribo_8$LogFC_RNA < -1 & datax_ribo_8$Sig == "Sig"] <- "DOWN"

############# SET separation ####################

genes_MITO<-datax_ribo_8$gene %in% result_df$Gene_Symbol
genes_SP<-datax_ribo_8$gene %in% list_genesSP$Symbol
genes_TM<-datax_ribo_8$gene %in% list_genesTM$Symbol

datax_ribo_8$Set <- NA
datax_ribo_8$Set[genes_MITO] <- "MITO"
datax_ribo_8$Set[genes_SP] <- "SP"
datax_ribo_8$Set[genes_TM] <- "TM"
categories_to_exclude <- c("MITO", "SP", "TM")
datax_ribo_8$Set <- ifelse(datax_ribo_8$Set %in% categories_to_exclude, datax_ribo_8$Set, "No SP/TM or mtRNA")
datax_ribo_8 <- datax_ribo_8[datax_ribo_8$Sig=="Sig",]

#combining TPM with significancy from Deseq for both conditions 4 and 8 hours


datax_ribo_4 <- data.frame(LogFC_RNA = ribo_results_4$log2FoldChange, Sig_RNA = as.numeric(ribo_results_4$padj), gene = ribo_results_4$gene)
datax_ribo_4$Sig_RNA[is.na(datax_ribo_4$Sig_RNA)] <- 1
temp_4 <- rep("NonSig", nrow(datax_ribo_4))
temp_4[datax_ribo_4$Sig_RNA < 0.05 & abs(datax_ribo_4$LogFC_RNA) > 1] <- "Sig"
datax_ribo_4$Sig <- factor(temp_4, levels = c("NonSig", "Sig"))

datax_ribo_4$diffexpressed <- "NO"
datax_ribo_4$diffexpressed[datax_ribo_4$LogFC_RNA < -1 & datax_ribo_4$Sig == "Sig"] <- "DOWN"
datax_ribo_4$diffexpressed[datax_ribo_4$LogFC_RNA > 1 & datax_ribo_4$Sig == "Sig"] <- "UP"

genes_MITO<-datax_ribo_4$gene %in% result_df$Gene_Symbol
genes_SP<-datax_ribo_4$gene %in% list_genesSP$Symbol
genes_TM<-datax_ribo_4$gene %in% list_genesTM$Symbol

datax_ribo_4$Set <- NA
datax_ribo_4$Set[genes_MITO] <- "MITO"
datax_ribo_4$Set[genes_SP] <- "SP"
datax_ribo_4$Set[genes_TM] <- "TM"
categories_to_exclude <- c("MITO", "SP", "TM")
datax_ribo_4$Set <- ifelse(datax_ribo_4$Set %in% categories_to_exclude, datax_ribo_4$Set, "No SP/TM or mtRNA")



####################
##################
##################
#############
##################
################
###############

datax_ribo_8_up <- datax_ribo_8[datax_ribo_8$diffexpressed=="UP",]
datax_ribo_8_down <- datax_ribo_8[datax_ribo_8$diffexpressed=="DOWN",]
common_genes_up <- datax_ribo_8_up[datax_ribo_8_up$gene %in% datax_ribo_4$gene,]
common_genes_down <- datax_ribo_8_down[datax_ribo_8_down$gene %in% datax_ribo_4$gene,]
datax_ribo_4$diffexpressed[datax_ribo_4$gene %in% common_genes_up$gene] <- "UP"
datax_ribo_4$diffexpressed[datax_ribo_4$gene %in% common_genes_down$gene] <- "DOWN"
common_genes <- datax_ribo_8[datax_ribo_8$gene %in% datax_ribo_4$gene,]
datax_ribo_4 <- datax_ribo_4[datax_ribo_4$gene %in% common_genes$gene,]


########## necessary from mitoanalysesRNA ############

################# MITOCHONDRIAL GENES ##################

# Read the list of genes with Entrez Gene IDs from the file
list_genesMITO <- read.table(file = "/Users/sgarcia/Downloads/mmilek-hdlbp_rev-c0e27b8/data/mitocarta2.txt", header = TRUE, sep = "\t")
# Define the Entrez Gene IDs you want to convert
ensembl_ids <- list_genesMITO$gene_id
# Load the Ensembl dataset for the human genome
ensembl <- useMart("ENSEMBL_MART_ENSEMBL", host = "www.ensembl.org", dataset = "hsapiens_gene_ensembl")
# Define your Ensembl dataset and attributes of interest
attributes_of_interest <- c("ensembl_gene_id", "external_gene_name")
# Create a list of unique Ensembl IDs
unique_ensembl_ids <- unique(ensembl_ids)
# Batch query to retrieve gene information for all unique Ensembl IDs
gene_info <- getBM(
  attributes = attributes_of_interest,
  filters = "ensembl_gene_id",
  values = unique_ensembl_ids,
  mart = ensembl
)

# Create a data frame to store the results
result_df <- data.frame(Ensembl_Gene_ID = unique_ensembl_ids, Gene_Symbol = NA)
# Fill in the Gene_Symbol column in result_df using the retrieved information
result_df$Gene_Symbol <- gene_info$external_gene_name[match(result_df$Ensembl_Gene_ID, gene_info$ensembl_gene_id)]

################# SP ####################
list_genesSP<-read.table(file="/Volumes/AG_Landthaler/sofia/riboSRP72/sp_list.txt", header = T, sep = "\t")
################# TM ####################
list_genesTM<-read.table(file="/Users/sgarcia/Downloads/mmilek-hdlbp_rev-c0e27b8/data/tm_list.txt", header = T, sep = "\t")
###################################### end ###############################

Rmembrane_genes_up_mine <- datax_ribo_8[datax_ribo_8$gene %in% membrane_genes & datax_ribo_8$diffexpressed=="UP" & datax_ribo_8$Sig=="Sig",]
Rcytosolic_genes_up_mine <- datax_ribo_8[datax_ribo_8$gene %in% cytosolic_genes & datax_ribo_8$diffexpressed=="UP" & datax_ribo_8$Sig=="Sig",]
Rundefined_genes_up_mine <-datax_ribo_8[datax_ribo_8$gene %in% undefined_genes & datax_ribo_8$diffexpressed=="UP" & datax_ribo_8$Sig=="Sig",]
RC_genes_paper_up <- datax_ribo_8[datax_ribo_8$gene %in% C_genes & datax_ribo_8$diffexpressed=="UP" & datax_ribo_8$Sig=="Sig",]
RM_genes_paper_up <- datax_ribo_8[datax_ribo_8$gene %in% E_genes & datax_ribo_8$diffexpressed=="UP" & datax_ribo_8$Sig=="Sig",]

nrow(Rmembrane_genes_up_mine)
nrow(Rcytosolic_genes_up_mine)
nrow(Rundefined_genes_up_mine)
nrow(RC_genes_paper_up)
nrow(RM_genes_paper_up)

###### Question: how many genes are cytosolic and membrane in my data and papers' data? ######

length(cytosolic_genes)
length(membrane_genes)
length(undefined_genes)
length(C_genes)
length(E_genes)

# Find the common elements between the two vectors
common_elements_cyto <- intersect(cytosolic_genes, C_genes)
common_elements_E <- intersect(membrane_genes, E_genes)
# Count the number of common elements
coincide_countC <- length(common_elements_cyto)
coincide_countE <- length(common_elements_E)
# Print the count of coinciding values
print(coincide_countC)
print(coincide_countE)

######### There is a upregulation of genes supposed to be cytosolic , 70 genes are significantly upregulated being cytosolic, only 8 veing ER########

####### next, how many of those genes have SP or TM? #######

############# 2 genes contain SP, 5 TM and 0 are mtRNA #######

SP_membrane_up <- list_genesSP[list_genesSP$Symbol %in% Rmembrane_genes_up_mine$gene,]
nrow(SP_membrane_up)
TM_membrane_up <- list_genesTM[list_genesTM$Symbol %in% Rmembrane_genes_up_mine$gene,]
nrow(TM_membrane_up)
MITO_membrane_up <- result_df[result_df$Gene_Symbol %in% Rmembrane_genes_up_mine$gene,]
nrow(MITO_membrane_up)

###### what are the cytosolic genes, there is only 1 with TM, the rest are 0 ########

SP_cytosolic_up <- list_genesSP[list_genesSP$Symbol %in% Rcytosolic_genes_up_mine$gene,]
nrow(SP_cytosolic_up)
TM_cytosolic_up <- list_genesTM[list_genesTM$Symbol %in% Rcytosolic_genes_up_mine$gene,]
nrow(TM_cytosolic_up)
MITO_cytosolic_up <- result_df[result_df$Gene_Symbol %in% Rcytosolic_genes_up_mine$gene,]
nrow(MITO_cytosolic_up)
undefined_up <- result_df[result_df$Gene_Symbol %in% Rundefined_genes_up_mine$gene,]
nrow(undefined_up)

#### 1. they dont contain SP or TM (TM only 1)

cytosolic_genes_TE<-gost(query = Rcytosolic_genes_up_mine$gene, organism = 'hsapiens',significant=T)
plot_cytosolic_TE<-gostplot(cytosolic_genes_TE,capped = TRUE, interactive = T)
cytosolic_TE_enrich<-enrichGO(Rcytosolic_genes_up_mine$gene, OrgDb='org.Hs.eg.db', keyType= "SYMBOL",ont= "MF",pAdjustMethod = "BH",pvalueCutoff  = 0.05,qvalueCutoff  = 0.05)
dotplot(cytosolic_TE_enrich)
barplot(cytosolic_TE_enrich, showCategory=10) 

datax_ribo4_up<-datax_ribo_4[datax_ribo_4$diffexpressed=="UP",]
datax_ribo8_down<-datax_ribo_8[datax_ribo_8$diffexpressed=="DOWN",]

membrane_TE_enrich<-enrichGO(datax_ribo4_up$gene, OrgDb='org.Hs.eg.db', keyType= "SYMBOL",ont= "MF",pAdjustMethod = "BH",pvalueCutoff  = 0.05,qvalueCutoff  = 0.05)
barplot(membrane_TE_enrich, showCategory=10)

membrane_TE_enrich<-enrichGO(datax_ribo4_down$gene, OrgDb='org.Hs.eg.db', keyType= "SYMBOL",ont= "ALL",pvalueCutoff  = 0.05,qvalueCutoff  = 0.05)
membrane_TE_enrich<-enrichGO(datax_ribo8_down$gene, OrgDb='org.Hs.eg.db', keyType= "SYMBOL",ont= "ALL",pvalueCutoff  = 0.05,qvalueCutoff  = 0.05)
barplot(membrane_TE_enrich,showCategory=10)

cytosolic_genes_TE<-gost(query = datax_ribo4_down$gene, organism = 'hsapiens',significant=T)

plot_cytosolic_TE<-gostplot(cytosolic_genes_TE,capped = TRUE, interactive = T)

datax_ribo_8_up$LogFC_RNA<-sort(datax_ribo_8_up$LogFC_RNA,decreasing=TRUE)

datax_ribo4_up_50<-datax_ribo_8_up[1:800,]

TE_enrich<-enrichGO(datax_ribo4_up_50$gene, OrgDb='org.Hs.eg.db', keyType= "SYMBOL",ont= "BP",pvalueCutoff  = 0.05,qvalueCutoff  = 0.05)
barplot(TE_enrich)

cnetplot(TE_enrich, node_label="all",color_category='firebrick', 
         color_gene='steelblue',cex_label_gene = 1,showCategory=8,categorySize="geneNum")

########## down ############

Rmembrane_genes_down_mine <- datax_ribo_8_down[datax_ribo_8_down$gene %in% membrane_genes,]

Rcytosolic_genes_down_mine <- datax_ribo_8_down[datax_ribo_8_down$gene %in% cytosolic_genes,]

cytosolic_TE_enrich_down<-enrichGO(Rcytosolic_genes_down_mine$gene, OrgDb='org.Hs.eg.db', keyType= "SYMBOL",ont= "CC",pAdjustMethod = "BH",pvalueCutoff  = 0.05,qvalueCutoff  = 0.05)
barplot(cytosolic_TE_enrich_down, showCategory=10) 

membrane_TE_enrich_down<-enrichGO(Rmembrane_genes_down_mine$gene, OrgDb='org.Hs.eg.db', keyType= "SYMBOL",ont= "BP",pAdjustMethod = "BH",pvalueCutoff  = 0.05,qvalueCutoff  = 0.05)
barplot(membrane_TE_enrich_down, showCategory=10) 


###############

datax_ribo_8$localisation <- NA
datax_ribo_8$localisation[datax_ribo_8$gene %in% membrane_genes] <- "membrane"
datax_ribo_8$localisation[datax_ribo_8$gene %in% cytosolic_genes] <- "cytosolic"
datax_ribo_8$localisation[datax_ribo_8$gene %in% undefined_genes] <- "undefined"


datax_ribo_4$localisation <- NA
datax_ribo_4$localisation[datax_ribo_4$gene %in% membrane_genes] <- "membrane"
datax_ribo_4$localisation[datax_ribo_4$gene %in% cytosolic_genes] <- "cytosolic"
datax_ribo_4$localisation[datax_ribo_4$gene %in% undefined_genes] <- "undefined"

genes_4_ribo <- datax_ribo_4$gene 
gene_info_4 <- getBM(
  attributes = c("external_gene_name", "gene_biotype"),
  filters = "external_gene_name",
  values = genes_4_ribo,
  mart = ensembl
)
datax_ribo_4 <- datax_ribo_4 %>%
  left_join(gene_info_4, by = c("gene" = "external_gene_name"))%>%
  mutate(gene_biotype = coalesce(gene_biotype, "noinfo"))

genes_8_ribo <- datax_ribo_8$gene 
gene_info_8 <- getBM(
  attributes = c("external_gene_name", "gene_biotype"),
  filters = "external_gene_name",
  values = genes_8_ribo,
  mart = ensembl
)
datax_ribo_8 <- datax_ribo_8 %>%
  left_join(gene_info_8, by = c("gene" = "external_gene_name"))%>%
  mutate(gene_biotype = coalesce(gene_biotype, "noinfo"))


# Add a column to datax_8 to indicate the source
datax_ribo_8 <- datax_ribo_8 %>%
  mutate(Source = "datax_ribo_8")

# Add a column to datax_4 to indicate the source
datax_ribo_4 <- datax_ribo_4 %>%
  mutate(Source = "datax_ribo_4")

# Combine the two data frames into one
combined_data_ribo <- bind_rows(datax_ribo_8, datax_ribo_4)
combined_data_ribo<- combined_data_ribo %>%
  mutate(Source = ifelse(Source == "datax_ribo_4", "4 hours", ifelse(Source == "datax_ribo_8", "8 hours", Source)))

# Assuming melted_data is your data frame
# Calculate the count of observations for each group

# Create the violin plot UP

combined_data_plot_ribo <-  combined_data_ribo
combined_data_plot_ribo$Group <- paste0(combined_data_plot_ribo$Set, ' (', combined_data_plot_ribo$Source, ')')
combined_data_plot_ribo <- combined_data_plot_ribo %>%
  filter(localisation %in% c("cytosolic", "membrane") & localisation != "undefined")
count_data_ribo <- combined_data_plot_ribo %>%
  group_by(Group,localisation) %>%
  summarise(Count = n(), .groups = 'drop')


ggplot(combined_data_plot_ribo, aes(x = Group, y = LogFC_RNA)) +
  geom_violin(aes(fill = Set),alpha=0.5) + geom_line(aes(group=gene), linewidth=0.4, color="gray") +
  geom_hline(yintercept = 0) + 
  geom_jitter(data=combined_data_plot_ribo, aes(x=Group, y=LogFC_RNA), width=.1, alpha=.5) +
  labs(title = "LogFC TE on ER", x = "Sets", y = "LogFC TE") +
  scale_fill_manual(values = c("SP" = "orange", "TM" = "dark green", "MITO" = "purple", "No SP/TM or mtRNA" = "grey")) +
  theme_minimal() +
  geom_text(data = count_data_ribo, aes(x = Group, y = Inf, label = Count), vjust = 1.5, size = 3, position = position_jitter(width = 0.1)) + 
  ylim(-5,7) + facet_grid(~localisation) + theme(
    axis.text.x = element_text(angle = 70, hjust = 1, size = 14), # Increase rotation angle and adjust text size
    axis.text.y = element_text(size = 14), # Increase axis text size
    axis.title.x = element_text(size = 20), # Increase axis title size
    axis.title.y = element_text(size = 20), # Increase axis title size
    plot.title = element_text(size = 20),   # Increase plot title size
    legend.text = element_text(size = 12),  # Increase legend text size
    legend.title = element_text(size = 16), # Increase legend title size
    strip.text = element_text(size = 16),   # Increase facet label size if needed
    panel.spacing = unit(2, "lines")
  )


####### simplified plot ##########

# Filter the data to include only "cytosolic" and "membrane" localizations
combined_data_ribo <- bind_rows(datax_ribo_8, datax_ribo_4)
combined_data_ribo<- combined_data_ribo %>%
  mutate(Source = ifelse(Source == "datax_ribo_4", "4 hours", ifelse(Source == "datax_ribo_8", "8 hours", Source)))

combined_data_ribo_plot <-  combined_data_ribo
filtered_data_ribo <- combined_data_ribo_plot %>%
  filter(localisation %in% c("cytosolic", "membrane") & localisation != "undefined")
count_data <- filtered_data %>%
  group_by(localisation,Source) %>%
  summarise(Count = n(), .groups = 'drop')


# Create the plot with the filtered data without using "Group"
ggplot(filtered_data_ribo, aes(x = Source, y = LogFC_RNA)) +
  geom_violin(aes(fill = localisation), alpha = 0.5) + 
  geom_line(aes(group = gene), linewidth = 0.4, color = "gray") +
  geom_hline(yintercept = 0) + 
  geom_jitter(data = filtered_data_ribo, aes(x = Source, y = LogFC_RNA), width = .1, alpha = .3) +
  labs(title = "LogFC TE on ER", x = "Time points", y = "LogFC TE") +
  scale_fill_manual(values = c("cytosolic" = "darkblue", "membrane" = "orange")) +
  theme_minimal() +
  geom_text(data = count_data, aes(x = Source, y = Inf, label = Count), vjust = 1.5, size = 3, position = position_jitter(width = 0.1)) + 
  ylim(-5, 5) +
  facet_wrap(~localisation, scales = "free") +  # Use facet_wrap to separate the plots by localization
  theme(
    axis.text.x = element_text(angle = 70, hjust = 1, size = 14), # Increase rotation angle and adjust text size
    axis.text.y = element_text(size = 14), # Increase axis text size
    axis.title.x = element_text(size = 20), # Increase axis title size
    axis.title.y = element_text(size = 20), # Increase axis title size
    plot.title = element_text(size = 20),   # Increase plot title size
    legend.text = element_text(size = 12),  # Increase legend text size
    legend.title = element_text(size = 16), # Increase legend title size
    strip.text = element_text(size = 16),   # Increase facet label size if needed
    panel.spacing = unit(2, "lines")
  )

############ get protein lengths #############

combined_data_cyto_protein <- combined_data_plot_ribo[combined_data_plot_ribo$localisation=="cytosolic" & combined_data_plot_ribo$Set=="No SP/TM or mtRNA" & combined_data_plot_ribo$Source=="8 hours" & combined_data_plot_ribo$diffexpressed=="UP",]
nrow(combined_data_cyto_protein)
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
gene_names <- combined_data_cyto_protein$gene
# Use the listAttributes function to get valid attribute names
valid_attributes <- listAttributes(ensembl)

gene_info <- getBM(attributes = c("external_gene_name", "cds_length","transcript_length","transcript_biotype"), filters = "external_gene_name", values = gene_names, mart = ensembl)
protein_coding_genes <- gene_info[gene_info$transcript_biotype=="protein_coding",]

protein_long_TE <- protein_coding_genes[protein_coding_genes$cds_length > 1000,]

# Extract unique gene names
unique_genes_TE <- unique(protein_long_TE$external_gene_name)

# Count the number of unique genes
unique_gene_count_TE <- length(unique_genes_TE)

### From 193 proteins, 180 of those are considered long proteins, as they have more than 500 aminoacids


library(writexl)

# Define the file path for your Excel file
file_path <- "/Users/sgarcia/Desktop/common_proteins.xlsx"

# Write the data frame to an Excel file
write_xlsx(protein_long_TE, file_path)

rnadata$gene <- rownames(rnadata)

common_proteins <- protein_long[protein_long$external_gene_name %in% protein_long_TE$external_gene_name,]
unique_proteins <- unique(common_proteins$external_gene_name)
length(unique_proteins)

unique_proteins_e<-enrichGO(unique_proteins, OrgDb='org.Hs.eg.db', keyType= "SYMBOL",ont= "BP",pAdjustMethod = "BH",pvalueCutoff  = 0.05,qvalueCutoff  = 0.05)

dotplot(unique_proteins_e)
barplot(unique_proteins_e, showCategory=10) 

# Write the data frame to an Excel file
write_xlsx(common_proteins, file_path)

unique(datax_ribo_8_up$gene)

####### lncRNAs #####

lncRNA_down<-datax_ribo_8[datax_ribo_8$diffexpressed=="DOWN" & datax_ribo_8$gene_biotype=="lncRNA",]
lncRNA_up<-datax_ribo_8[datax_ribo_8$diffexpressed=="UP" & datax_ribo_8$gene_biotype=="lncRNA",]

ribo_up<-enrichGO(datax_ribo_8_up$gen, OrgDb='org.Hs.eg.db', keyType= "SYMBOL",ont= "MF",pAdjustMethod = "BH",pvalueCutoff  = 0.05,qvalueCutoff  = 0.05)

dotplot(ribo_up)

ribo_down<-enrichGO(datax_ribo_8_down$gene, OrgDb='org.Hs.eg.db', keyType= "SYMBOL",ont= "CC",pAdjustMethod = "BH",pvalueCutoff  = 0.05,qvalueCutoff  = 0.05)

dotplot(ribo_down)

down_down<-datax_8_down[list_genesSP$Symbol %in% datax_ribo_8_down$gene,]



############## heatmap ##################

library(pheatmap)
library(ComplexHeatmap)
dataxtry_ribo<-as.data.frame(cbind(datax_ribo_4,datax_ribo_8))
dataxtry_RNA<-as.data.frame(cbind(datax_4,datax_8))

# Sort the data based on LogFC_RNA values from datax_ribo_8 (10th column)
sorted_dataxtry_ribo <- dataxtry_ribo[order(dataxtry_ribo[,7], decreasing = F), ]#dependiendo en si sumas las tablas mas adelante o no, la seleccion de la columna varia entre 5 o 10
sorted_dataxtry_RNA <- dataxtry_RNA[order(dataxtry_RNA[,6], decreasing = F), ]

# Merge the two data frames using a full outer join and fill missing values with 0
y <- merge(sorted_dataxtry_ribo, sorted_dataxtry_RNA, by = 3, all = TRUE, fill = 0)

# Replace NA values with 0 if necessary
y[is.na(y)] <- 0


# Subset data for the top 20 genes for heatmap
heatmap_data_20 <- sorted_dataxtry_ribo[1:20, c(1, 7)]
colnames(heatmap_data_20) <- c("4 hours", "8 hours")

# Set rownames to gene symbols for top 20 genes
rownames(heatmap_data_20) <- sorted_dataxtry$gene[1:20]

# Create annotation for top 20 genes
annotation_row_20 <- data.frame(Localisation = sorted_dataxtry$localisation[1:20])
annotation_row_20$Localisation <- factor(annotation_row_20$Localisation)

# Define the colors for the heatmap and for the annotations
heatmap_colors <- colorRamp2(c(-3, 0, 3), c("blue", "white", "red"))
ann_colors <- list(Localisation = c("cytosolic" = "red", "membrane" = "blue", "undefined" = "yellow", "MITO" = "green"))

# Row annotations
ra = rowAnnotation(df = annotation_row_20, col = ann_colors)

# Generate the heatmap
Heatmap(heatmap_data_20, name = "LogFC", 
        col = heatmap_colors, 
        row_names_gp = gpar(fontsize = 10),
        row_title = "Genes", 
        column_title = "Time Points",
        show_row_names = TRUE, 
        show_column_names = TRUE,
        heatmap_legend_param = list(title = "Log Fold Change"),
        row_annotation = ra)

######## large plot #########

### ribo

heatmap_data_all_ribo <- dataxtry_ribo[, c(1, 7)]
colnames(heatmap_data_all_ribo) <- c("4 hours Ribo", "8 hours Ribo")

# Adjusting margin to make space for colorbar title
colorbar_margin <- c(0, 4, 0, 4)  # top, right, bottom, left
heatmap_data_all_ribo<-as.matrix(heatmap_data_all_ribo)
# Generate the heatmap for all genes without annotation column
p7<-pheatmap::pheatmap(heatmap_data_all_ribo, 
         show_rownames = FALSE,
         color = colorRampPalette(c("lightblue", "black", "yellow"))(200),
         cluster_cols = TRUE,
         colorbar_title = "Log Fold Change",
         main = "Significantly changing genes in ER Ribo-Seq data",  # This is optional; just if you want a title for the entire heatmap
         mar = colorbar_margin)

#### rna

heatmap_data_all_RNA <- dataxtry_RNA[, c(1, 6)]
colnames(heatmap_data_all_RNA) <- c("4 hours RNA","8 hours RNA")

# Adjusting margin to make space for colorbar title
colorbar_margin <- c(0, 4, 0, 4)  # top, right, bottom, left
heatmap_data_all_RNA<-as.matrix(heatmap_data_all_RNA)
# Generate the heatmap for all genes without annotation column
p8<-pheatmap::pheatmap(heatmap_data_all_RNA, 
                   show_rownames = FALSE,
                   color = colorRampPalette(c("lightblue", "black", "yellow"))(200),
                   cluster_cols = TRUE,
                   colorbar_title = "Log Fold Change",
                   main = "Significantly changing genes in ER RNA-Seq data",  # This is optional; just if you want a title for the entire heatmap
                   mar = colorbar_margin)

