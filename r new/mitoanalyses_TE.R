library(DESeq2)
library(EnsDb.Hsapiens.v86)
library(dplyr)
library(biomaRt)
library(ggplot2)

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


# path
dir.path <- "/Volumes/AG_Landthaler/amanukyan/Projects/RiboSeq/Pipeline/run/sofia/output/FeatureCounts"
getwd()
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
temp <- rep("NonSig", nrow(datax_ribo_8))
temp[datax_ribo_8$Sig_RNA < 0.05 & abs(datax_ribo_8$LogFC_RNA) > 1] <- "Sig"
datax_ribo_8$Sig <- factor(temp, levels = c("NonSig", "Sig"))

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

datax_ribo_4$Set <- NA
datax_ribo_4$Set[genes_MITO] <- "MITO"
datax_ribo_4$Set[genes_SP] <- "SP"
datax_ribo_4$Set[genes_TM] <- "TM"
categories_to_exclude <- c("MITO", "SP", "TM")
datax_ribo_4$Set <- ifelse(datax_ribo_4$Set %in% categories_to_exclude, datax_ribo_4$Set, "No SP/TM or mtRNA")


#combining TPM with significancy from Deseq for both conditions 4 and 8 hours
datax_ribo_4 <- data.frame(LogFC_RNA = ribo_results_4$log2FoldChange, Sig_RNA = as.numeric(ribo_results_4$padj), gene = ribo_results_4$gene)
datax_ribo_4$Sig_RNA[is.na(datax_ribo_4$Sig_RNA)] <- 1
temp <- rep("NonSig", nrow(datax_ribo_4))
temp[datax_ribo_4$Sig_RNA < 0.05 & abs(datax_ribo_4$LogFC_RNA) > 1] <- "Sig"
datax_ribo_4$Sig <- factor(temp, levels = c("NonSig", "Sig"))

datax_ribo_4$diffexpressed <- "NO"
datax_ribo_4$diffexpressed[datax_ribo_4$LogFC_RNA < -1 & datax_ribo_4$Sig == "Sig"] <- "DOWN"
genes_MITO<-datax_ribo_4$gene %in% result_df$Gene_Symbol
genes_SP<-datax_ribo_4$gene %in% list_genesSP$Symbol
genes_TM<-datax_ribo_4$gene %in% list_genesTM$Symbol


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


# Add a column to datax_8 to indicate the source
datax_ribo_8 <- datax_ribo_8 %>%
  mutate(Source = "datax_ribo_8")
# Add a column to datax_4 to indicate the source
datax_ribo_4 <- datax_ribo_4 %>%
  mutate(Source = "datax_ribo_4")

# Combine the two data frames into one
combined_data_ribo <- bind_rows(datax_ribo_8, datax_ribo_4)
# Create the violin plot UP
combined_data_up <- combined_data_ribo[combined_data_ribo$diffexpressed == "UP",]
combined_data_up_plot <-  combined_data_up
combined_data_up_plot$Group <- paste0(combined_data_up_plot$Set, ' (', combined_data_up_plot$Source, ')')
count_data_up <- combined_data_up_plot %>%
  group_by(Group) %>%
  summarise(Count = n(), .groups = 'drop')



ggplot(combined_data_up_plot, aes(x = Group, y = LogFC_RNA)) +
  geom_violin(aes(fill = Set)) + geom_line(aes(group=gene), linewidth=0.4, color="gray") +
  geom_hline(yintercept = 0) + 
  geom_jitter(data=combined_data_up_plot, aes(x=Group, y=LogFC_RNA), width=.1, alpha=.5) +
  labs(title = "LogFC RNA Downregulated on ER", x = "Gene sets", y = "LogFC RNA") +
  scale_fill_manual(values = c("SP" = "blue", "TM" = "green", "MITO" = "red", "Cytosol" = "grey")) +
  theme_minimal() +
  geom_text(data = count_data_up, aes(x = Group, y = Inf, label = Count), vjust = 1.5, size = 3, position = position_jitter(width = 0.1)) + 
  ylim(-3,7)

# Create the violin plot DOWN

combined_data_down <- combined_data_ribo[combined_data_ribo$diffexpressed == "DOWN",]

combined_data_down_plot <-  combined_data_down
combined_data_down_plot$Group <- paste0(combined_data_down_plot$Set, ' (', combined_data_down_plot$Source, ')')

count_data_down <- combined_data_down_plot %>%
  group_by(Group) %>%
  summarise(Count = n(), .groups = 'drop')

ggplot(combined_data_down_plot, aes(x = Group, y = LogFC_RNA)) +
  geom_violin(aes(fill = Set)) + geom_line(aes(group=gene), linewidth=0.4, color="gray") +
  geom_hline(yintercept = 0) + 
  geom_jitter(data=combined_data_down_plot, aes(x=Group, y=LogFC_RNA), width=.1, alpha=.5) +
  labs(title = "LogFC RNA Downregulated on ER", x = "Gene sets", y = "LogFC RNA") +
  scale_fill_manual(values = c("SP" = "blue", "TM" = "green", "MITO" = "red", "Cytosol" = "grey")) +
  theme_minimal() +
  geom_text(data = count_data_down, aes(x = Group, y = Inf, label = Count), vjust = 1.5, size = 3, position = position_jitter(width = 0.1)) + 
  ylim(-7,1)

