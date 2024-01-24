library(DESeq2)

#needed for riborex
#  install.packages('devtools')
#library(devtools)
#options(unzip='internal')
#devtools::install_github('smithlabcode/riborex')

library(edgeR)
library(fdrtool)
library(riborex)
library(EnsDb.Hsapiens.v86)
library(dplyr)
library(biomaRt)
library(ggplot2)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(clusterProfiler)
#(gage)
library(gprofiler2)
library(knitr)
library(enrichplot)
library(DOSE) # needed to convert to enrichResult object
library(gridExtra)
library(ggrepel)
library(reshape2)
library(tidyverse)

# library(riborex)
#source("/Volumes/sofia/DESeq2Rex.R")

####
# Import data ####
####

## replicate E0R2 from RNA is excluded as seems outlier in PCA, corrplots and has lower counts ####
#### Replicate E81 was also excluded as has contamination of UTR and introns and corrplot looked not so good ####

# path
dir.path <- "/Volumes/AG_Landthaler/amanukyan/Projects/RiboSeq/Pipeline/run/sofia/output/FeatureCounts"
getwd()
# human data 
list_files <- list.files(dir.path)
list_files <- list_files[!grepl("summary", list_files)]

rnadata <- read.table(file="/Volumes/AG_Landthaler/amanukyan/Projects/RiboSeq/Pipeline/run/sofia_rna/output/count.txt", header = T, sep = "\t")
colnames(rnadata) <- lapply(colnames(rnadata), function(x) strsplit(x, split = "_S[0-9]+_")[[1]][1])
rownames(rnadata) <- rnadata[,1]
rnadata <- rnadata[,-1]
rnadata <- rnadata[,c(14,16,17:19,20:21)]
colnames(rnadata) <- gsub("output.bam_files.", "", colnames(rnadata))

####
# DE analysis ####
####

# remove low count genes
count_gene <- apply(rnadata, 1, max) 
rnadata <- rnadata[which(count_gene > 10),]

# control treatment
rnaCond <- c("ctrl","ctrl","time 4 hours","time 4 hours","time 4 hours","time 8 hours","time 8 hours")
rnaCond<-factor(rnaCond)
# get riborex data
print(head(rnadata))

#Differential expression of RNA on ER among conditions
dds <- DESeqDataSetFromMatrix(countData = rnadata,
                              DataFrame(rnaCond),
                              design= ~ rnaCond)
dds <- DESeq(dds)
unique_rnaCond<- as.character(rev(unique(rnaCond)))
rna_results <- NULL
for(i in 1:(length(unique_rnaCond)-1)){
  for(j in (i+1):(length(unique_rnaCond))){
    comparison <- c("rnaCond", unique_rnaCond[i], unique_rnaCond[j])
    print(comparison)
    cur_results <- as.data.frame(results(dds, contrast = comparison))
    rna_results <- rbind(rna_results, data.frame(cur_results, gene = rownames(cur_results), comparison = paste0(comparison[2:3], collapse = "_vs_")))
  }
}
rna_results$padj[is.na(rna_results$padj)] <- 1
rna_results_8<-filter(rna_results, comparison == "time 8 hours_vs_ctrl")
rna_results_4<-filter(rna_results, comparison == "time 4 hours_vs_ctrl")

# Check for matching row names
matching_rows <- rownames(rnadata) %in% rna_results_8$gene

# Print the result
print(matching_rows)

# Select rows from table1 where matching_rows is TRUE
selected_rows <- rnadata[matching_rows, ]

# Print the selected rows
print(selected_rows)

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

####
# Visualize ####
####

######RNA cond2 vs cond1, 

#calculation of the TPM log mean will be used to see which genes are up or downregulated. DDSeq analises will be useful to determine which of those
#are significant and which not. But the numeric value we get ourselves performing our own TPM normalization.

datax_4 <- data.frame(LogFC_RNA = rna_results_4$log2FoldChange, Sig_RNA = as.numeric(rna_results_4$padj), gene = rna_results_4$gene)
datax_4$Sig_RNA[is.na(datax_4$Sig_RNA)] <- 1
temp <- rep("NonSig", nrow(datax_4))
temp[datax_4$Sig_RNA < 0.05 & abs(datax_4$LogFC_RNA) > 1] <- "Sig"
datax_4$Sig <- factor(temp, levels = c("NonSig", "Sig"))
cur_LogMean_RNA_Data_4 <- LogMean_RNA_Data[datax_4$gene,]


datax_8 <- data.frame(LogFC_RNA = rna_results_8$log2FoldChange, Sig_RNA = as.numeric(rna_results_8$padj), gene = rna_results_8$gene)
datax_8$Sig_RNA[is.na(datax_8$Sig_RNA)] <- 1
temp <- rep("NonSig", nrow(datax_8))
temp[datax_8$Sig_RNA < 0.05 & abs(datax_8$LogFC_RNA) > 1] <- "Sig"
datax_8$Sig <- factor(temp, levels = c("NonSig", "Sig"))
cur_LogMean_RNA_Data_8 <- LogMean_RNA_Data[datax_8$gene,]

# Check for NA values in table1
rows_with_na <- rowSums(is.na(cur_LogMean_RNA_Data_8)) > 0

# Select rows from table1 where row names have NA values
selected_rows <- cur_LogMean_RNA_Data_8[rows_with_na, ]

# Print the selected rows
print(selected_rows)

#condition 8 hours is column 3, condition 4 hours is column 2, ctrl is column 1
#x and y are normalized data

cur_LogMean_RNA_Data_4_select <- data.frame(x = as.numeric(cur_LogMean_RNA_Data_4[,1]), 
                                     y = as.numeric(cur_LogMean_RNA_Data_4[,2]), 
                                     Sig = as.character(datax_4$Sig))


cur_LogMean_RNA_Data_8_select <- data.frame(x = as.numeric(cur_LogMean_RNA_Data_8[,1]), 
                                     y = as.numeric(cur_LogMean_RNA_Data_8[,3]), 
                                     Sig = as.character(datax_8$Sig))

#### Addition of extra columns and information about significance, mixing our own TPM calculation and Desq calculations.

#### For the condition of 4 hours

# add a column of NAs
cur_LogMean_RNA_Data_4_select$diffexpressed <- "NO"
# if log2Foldchange > 1.5 and pvalue < 0.05, set as "UP" 
cur_LogMean_RNA_Data_4_select$diffexpressed[datax_4$LogFC_RNA > 1 & datax_4$Sig == "Sig"] <- "UP"
# if log2Foldchange < -1 and pvalue < 0.05, set as "DOWN"
cur_LogMean_RNA_Data_4_select$diffexpressed[datax_4$LogFC_RNA < -1 & datax_4$Sig == "Sig"] <- "DOWN"
cur_LogMean_RNA_Data_4_select$gene_symbol<-datax_4$gene
cur_LogMean_RNA_Data_4_select$delabel <- NA
cur_LogMean_RNA_Data_4_select$delabel[cur_LogMean_RNA_Data_4_select$diffexpressed != "NO"] <- cur_LogMean_RNA_Data_4_select$gene_symbol[cur_LogMean_RNA_Data_4_select$diffexpressed != "NO"]

### For the condition of 8 hours

# add a column of NAs
cur_LogMean_RNA_Data_8_select$diffexpressed <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
cur_LogMean_RNA_Data_8_select$diffexpressed[datax_8$LogFC_RNA > 1 & datax_8$Sig == "Sig"] <- "UP"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
cur_LogMean_RNA_Data_8_select$diffexpressed[datax_8$LogFC_RNA < -1 & datax_8$Sig == "Sig"] <- "DOWN"
cur_LogMean_RNA_Data_8_select$gene_symbol<-datax_8$gene
cur_LogMean_RNA_Data_8_select$delabel <- NA
cur_LogMean_RNA_Data_8_select$delabel[cur_LogMean_RNA_Data_8_select$diffexpressed != "NO"] <- cur_LogMean_RNA_Data_8_select$gene_symbol[cur_LogMean_RNA_Data_8_select$diffexpressed != "NO"]

##### Selection here of genes that are significantly downregulated from the mother table


genes_down_RNA_4<-cur_LogMean_RNA_Data_4_select%>%
  dplyr::filter(diffexpressed == "DOWN")%>%
  dplyr::filter(Sig=="Sig")%>%
  pull(gene_symbol)

genes_down_RNA_8<-cur_LogMean_RNA_Data_8_select%>%
  dplyr::filter(diffexpressed == "DOWN")%>%
  dplyr::filter(Sig=="Sig")%>%
  pull(gene_symbol)

####### SP vs downregulated ##########

#I want to test whether the downregulation of SP genes is progressive by looking at the table of time 4 and 8 and seeing 
#list of genes with signal peptide

list_genesER<-read.table(file="/Volumes/AG_Landthaler/sofia/riboSRP72/sp_list.txt", header = T, sep = "\t")

#next are the genes with SP present on the downregulated genes in the RNA data after 8 h depletion, they are 58 genes same after 4 and 8 hours

genes_SP_down_8<-list_genesER[list_genesER$Symbol %in% genes_down_RNA_8,]

genes_SP_down_8<-print(genes_SP_down_8$Symbol)

#when we dont apply sig up in the filtering of the downregulated genes, 58 genes containing SP from 339 are included in this group of downregulated genes
#it is reduced to 13 when applying Sig.


sp<-enrichGO(genes_SP_down_8, OrgDb='org.Hs.eg.db', keyType= "SYMBOL",ont= "ALL",pAdjustMethod = "BH",pvalueCutoff  = 0.05,qvalueCutoff  = 0.05)

dotplot(sp)

ue<-gost(query = genes_SP_down_8, organism = 'hsapiens',significant=T)

ue2<-gostplot(ue,capped = TRUE, interactive = F)

publish_gostplot(ue2,highlight_terms =c("GO:0016020","GO:0012505","GO:0098827","GO:0070062","GO:0005783","GO:0005576","GO:0007115","GO:0004714","KEGG:04514","REAC:R-HSA-1474244","MF:M07289"))

#####################
####################
####################


###### SP genes enrichment ######

question<-enrichGO(list_genesER$Symbol, OrgDb='org.Hs.eg.db', keyType= "SYMBOL",ont= "ALL",pAdjustMethod = "BH",pvalueCutoff  = 0.05,qvalueCutoff  = 0.05)

gseDO(list_genesER$Symbol)

dotplot(question)

#same but with 4 hours


genes_SP_down_4<-list_genesER[list_genesER$Symbol %in% genes_down_RNA_4,]

genes_SP_down_4<-print(genes_SP_down_4$Symbol)



##### TM list #######

list_genesTM<-read.table(file="/Users/sgarcia/Downloads/mmilek-hdlbp_rev-c0e27b8/data/tm_list.txt", header = T, sep = "\t")

genes_TM_down_8<-list_genesTM[list_genesTM$Symbol %in% genes_down_RNA_8,]

genes_TM_down_8<-print(genes_TM_down_8$Symbol)


tm_down_8<-enrichGO(genes_TM_down_8, OrgDb='org.Hs.eg.db', keyType= "SYMBOL",ont= "ALL",pAdjustMethod = "BH",pvalueCutoff  = 0.05,qvalueCutoff  = 0.05)

dotplot(tm_down_8)

tm_down_8_gost<-gost(query = genes_TM_down_8, organism = 'hsapiens',significant=T)

gostplot_TM_down_8<-gostplot(tm_down_8_gost,capped = TRUE, interactive = F)

publish_gostplot(gostplot_TM_down_8,highlight_terms =c("GO:0016020","GO:0012505","GO:0098827","GO:0070062","GO:0005783","GO:0005576","GO:0007115","GO:0004714","GO:0019199","GO:0005539","GO:0007155","GO:0070062","GO:0005783","TF:M02089","TF:M10115","TF:M08867","TF:M03557"))



########## TOTAL FRACTION #################

###### I eliminated Ctrl 3 as in PCA seems a bit outlier 

total_fractions<- read.table(file="/Volumes/AG_Landthaler/amanukyan/Projects/RiboSeq/Pipeline/run/sofia_rna/output/count.txt", header = T, sep = "\t")
colnames(total_fractions) <- lapply(colnames(total_fractions), function(x) strsplit(x, split = "_S[0-9]+_")[[1]][1])
rownames(total_fractions) <- total_fractions[,1]
total_fractions <- total_fractions[,c(23,24,26:30)]
colnames(total_fractions) <- gsub("output.bam_files.", "", colnames(total_fractions))
colnames(total_fractions) <- gsub(".bam", "", colnames(total_fractions))


####
# DE analysis ####
####

# remove low count genes
count_gene_total <- apply(total_fractions, 1, max) 
rnadata_total <- total_fractions[which(count_gene_total > 10),]

rnaCond_total <- c("ctrl","ctrl","time 4 hours","time 4 hours","time 4 hours","time 8 hours","time 8 hours")
rnaCond_total<-factor(rnaCond_total)
# get riborex data
print(head(rnadata_total))

# get regular deseq2 results
dds <- DESeqDataSetFromMatrix(countData = rnadata_total,
                              DataFrame(rnaCond_total),
                              design= ~ rnaCond_total)
dds <- DESeq(dds)
unique_rnaCond_total <- as.character(rev(unique(rnaCond_total)))
rna_results_total <- NULL
for(i in 1:(length(unique_rnaCond_total)-1)){
  for(j in (i+1):(length(unique_rnaCond_total))){
    comparison <- c("rnaCond_total", unique_rnaCond_total[i], unique_rnaCond_total[j])
    print(comparison)
    cur_results <- as.data.frame(results(dds, contrast = comparison))
    rna_results_total <- rbind(rna_results_total, data.frame(cur_results, gene = rownames(cur_results), comparison = paste0(comparison[2:3], collapse = "_vs_")))
  }
}
rna_results_total$padj[is.na(rna_results_total$padj)] <- 1
rna_results_total_8<-filter(rna_results_total, comparison == "time 8 hours_vs_ctrl")
rna_results_total_4<-filter(rna_results_total, comparison == "time 4 hours_vs_ctrl")

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
len_rna_total <- len[match(rownames(rnadata_total), len$hgnc_symbol),]
len_rna_total <- na.omit(len_rna_total)
rnadata_tpm_total <- rnadata_total[intersect(rownames(rnadata_total), len_rna_total$hgnc_symbol),]
rnadata_tpm_total <- sweep(rnadata_tpm_total * 1000, 1, len_rna_total$len, FUN = "/")
rnadata_tpm_total <- sweep(rnadata_tpm_total * 1000000, 2, colSums(rnadata_tpm_total), FUN = "/")

# calculate log RNA mean for each condition
cond_names <- rnaCond_total
Mean_RNA_Data_total <- aggregate(t(rnadata_tpm_total), list(cond_names), mean)
Mean_RNA_Data_total <- t(Mean_RNA_Data_total)
colnames(Mean_RNA_Data_total) <- Mean_RNA_Data_total[1,]
Mean_RNA_Data_total <- data.frame(apply(Mean_RNA_Data_total[-1,],2,as.numeric), row.names = rownames(Mean_RNA_Data_total[-1,]))
LogMean_RNA_Data_total <- log(Mean_RNA_Data_total + 0.001)

datax_total_4 <- data.frame(LogFC_RNA = rna_results_total_4$log2FoldChange, Sig_RNA = as.numeric(rna_results_total_4$padj), gene = rna_results_total_4$gene)
datax_total_4$Sig_RNA[is.na(datax_total_4$Sig_RNA)] <- 1
temp_total <- rep("NonSig", nrow(datax_total_4))
temp_total[datax_total_4$Sig_RNA < 0.05 & abs(datax_total_4$LogFC_RNA) > 1] <- "Sig"
datax_total_4$Sig <- factor(temp_total, levels = c("NonSig", "Sig"))
cur_LogMean_RNA_Data_total_4 <- LogMean_RNA_Data_total[datax_total_4$gene,]

cur_LogMean_RNA_Data_total_4_select <- data.frame(x = as.numeric(cur_LogMean_RNA_Data_total_4[,1]), 
                                                  y = as.numeric(cur_LogMean_RNA_Data_total_4[,3]), 
                                                  Sig = as.character(datax_total_4$Sig))

datax_total_8 <- data.frame(LogFC_RNA = rna_results_total_8$log2FoldChange, Sig_RNA = as.numeric(rna_results_total_8$padj), gene = rna_results_total_8$gene)
datax_total_8$Sig_RNA[is.na(datax_total_8$Sig_RNA)] <- 1
temp_total <- rep("NonSig", nrow(datax_total_8))
temp_total[datax_total_8$Sig_RNA < 0.05 & abs(datax_total_8$LogFC_RNA) > 1] <- "Sig"
datax_total_8$Sig <- factor(temp_total, levels = c("NonSig", "Sig"))
cur_LogMean_RNA_Data_total_8 <- LogMean_RNA_Data_total[datax_total_8$gene,]

cur_LogMean_RNA_Data_total_8_select <- data.frame(x = as.numeric(cur_LogMean_RNA_Data_total_8[,1]), 
                                           y = as.numeric(cur_LogMean_RNA_Data_total_8[,3]), 
                                           Sig = as.character(datax_total_8$Sig))

# add a column of NAs
cur_LogMean_RNA_Data_total_4_select$diffexpressed <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
cur_LogMean_RNA_Data_total_4_select$diffexpressed[datax_total_4$LogFC_RNA > 1 & datax_total_4$Sig == "Sig"] <- "UP"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
cur_LogMean_RNA_Data_total_4_select$diffexpressed[datax_total_4$LogFC_RNA < -1 & datax_total_4$Sig == "Sig"] <- "DOWN"

cur_LogMean_RNA_Data_total_4_select$gene_symbol<-datax_total_4$gene
cur_LogMean_RNA_Data_total_4_select$delabel <- NA
cur_LogMean_RNA_Data_total_4_select$delabel[cur_LogMean_RNA_Data_total_4_select$diffexpressed != "NO"] <- cur_LogMean_RNA_Data_total_4_select$gene_symbol[cur_LogMean_RNA_Data_total_4_select$diffexpressed != "NO"]

genes_down_RNA_total_4<-cur_LogMean_RNA_Data_total_4_select%>%
  dplyr::filter(diffexpressed == "DOWN")%>%
  dplyr::filter(Sig=="Sig")%>%
  select(gene_symbol)

# add a column of NAs
cur_LogMean_RNA_Data_total_8_select$diffexpressed <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
cur_LogMean_RNA_Data_total_8_select$diffexpressed[datax_total_8$LogFC_RNA > 1 & datax_total_8$Sig == "Sig"] <- "UP"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
cur_LogMean_RNA_Data_total_8_select$diffexpressed[datax_total_8$LogFC_RNA < -1 & datax_total_8$Sig == "Sig"] <- "DOWN"

cur_LogMean_RNA_Data_total_8_select$gene_symbol<-datax_total_8$gene
cur_LogMean_RNA_Data_total_8_select$delabel <- NA
cur_LogMean_RNA_Data_total_8_select$delabel[cur_LogMean_RNA_Data_total_8_select$diffexpressed != "NO"] <- cur_LogMean_RNA_Data_total_8_select$gene_symbol[cur_LogMean_RNA_Data_total_8_select$diffexpressed != "NO"]

genes_down_RNA_total_8<-cur_LogMean_RNA_Data_total_8_select%>%
  dplyr::filter(diffexpressed == "DOWN")%>%
  dplyr::filter(Sig=="Sig")%>%
  select(gene_symbol)

### visualization of the enrichment of the 23 genes downregulated on the total fraction
genes_down_8<-enrichGO(genes_down_RNA_total_8$gene_symbol, OrgDb='org.Hs.eg.db', keyType= "SYMBOL",ont= "CC",pAdjustMethod = "BH",pvalueCutoff  = 0.05,qvalueCutoff  = 0.05)
barplot(genes_down_8)

total_down_8_gost<-gost(query = genes_down_RNA_total_8$gene_symbol, organism = 'hsapiens',significant=T)
plot_total_gost<-gostplot(total_down_8_gost,capped = TRUE, interactive = T)

### up
genes_up_RNA_total_4<-cur_LogMean_RNA_Data_total_4_select%>%
  dplyr::filter(diffexpressed == "UP")%>%
  dplyr::filter(Sig=="Sig")%>%
  select(gene_symbol)

genes_up_RNA_total_8<-cur_LogMean_RNA_Data_total_8_select%>%
  dplyr::filter(diffexpressed == "UP")%>%
  dplyr::filter(Sig=="Sig")%>%
  select(gene_symbol)

genes_up_RNA_total_4[genes_up_RNA_total_4$gene_symbol %in% genes_up_RNA_total_8$gene_symbol,]

genes_up_8<-enrichGO(genes_up_RNA_total_8$gene_symbol, OrgDb='org.Hs.eg.db', keyType= "SYMBOL",ont= "BP",pAdjustMethod = "BH",pvalueCutoff  = 0.05,qvalueCutoff  = 0.05)
barplot(genes_up_8)

genes_up_8<-as.data.frame(genes_up_8)
# Extracting the corresponding GO terms from the ID column
selected_GO_terms <- genes_up_8$ID[1:5]
# Printing the selected GO terms
print(selected_GO_terms)

list_genes<-list(genes_up_RNA_total_4,genes_up_RNA_total_8)
total_up_gost<-gost(query = c(list_genes, query_1="4 hours treatment",query_2="8 hours treatment"),organism = 'hsapiens',significant=T,highlight=TRUE)
plot_total_gost<-gostplot(total_up_gost,interactive=F,capped=TRUE)
publish_gostplot(plot_total_gost, highlight_terms= c("GO:0006457","GO:0061077" ,"GO:0006986" ,"GO:0035966" ,"GO:0009408","TF:M11655_1","TF:M07250_1","GO:1903561","GO:0060590",
                                                     "GO:0005925","GO:0044183","GO:0031072","REAC:R-HSA-5683057","GO:0005925","MIRNA:hsa-miR-16-5p"))

#list of genes with signal peptide (339 genes)

list_genesER<-read.table(file="/Volumes/AG_Landthaler/sofia/riboSRP72/sp_list.txt", header = T, sep = "\t")

genes_SP_down_total<-list_genesER[list_genesER$Symbol %in% genes_down_RNA_total_8,]

genes_SP_down_total<-print(genes_SP_down_total$Symbol)

sp_total<-enrichGO(genes_SP_down_total, OrgDb='org.Hs.eg.db', keyType= "SYMBOL",ont= "ALL",pAdjustMethod = "BH",pvalueCutoff  = 0.05,qvalueCutoff  = 0.05)

dotplot(sp_total)

total_down_SP<-gost(query = genes_SP_down_total, organism = 'hsapiens',significant=T)

total_enrich_SP<-gostplot(total_down_SP,capped = TRUE, interactive = T)

####these genes SP depleted on the total fraction do not show any enrichment

###### TM list #####

list_genesTM<-read.table(file="/Users/sgarcia/Downloads/mmilek-hdlbp_rev-c0e27b8/data/tm_list.txt", header = T, sep = "\t")

genes_TM_down_8_total<-list_genesTM[list_genesTM$Symbol %in% genes_down_RNA_total_8,]

genes_TM_down_8_total<-print(genes_TM_down_8_total$Symbol)

#what we observe is that 297 from 339 genes with SP are downregulated in the total farction contrary to the ER fraction,
#when we apply sig in the filtering a bit up, the number goes down to 44.


####### comparison between the genes downregulated in ER and in total ########

comparison_ER_total_8 <- genes_down_RNA_total_8[genes_down_RNA_total_8$gene_symbol %in% datax_8_down$gene,]

comparison_ER_total_4<- genes_down_RNA_total_4[genes_down_RNA_total_4$gene_symbol %in% datax_4_down$gene,]

sis <- enrichGO(comparison_ER_total, OrgDb='org.Hs.eg.db', keyType= "SYMBOL",ont= "ALL",pAdjustMethod = "BH",pvalueCutoff  = 0.05,qvalueCutoff  = 0.05)

dotplot(sis)

non_comparison_ER_total <- genes_down_RNA_total_8[!genes_down_RNA_total_8 %in% genes_down_RNA_8]

non <- enrichGO(non_comparison_ER_total, OrgDb='org.Hs.eg.db', keyType= "SYMBOL",ont= "ALL",pAdjustMethod = "BH",pvalueCutoff  = 0.05,qvalueCutoff  = 0.05)

dotplot(non)

sis<-gost(query = non_comparison_ER_total, organism = 'hsapiens',significant=T)

sis2<-gostplot(sis,capped = TRUE, interactive = T)



#### From the comparison we observe how in ER fraction after 8 hours there is way more significant downregulation of genes than in the total fraction
#### What we can think of is that mRNAs are mislocalized to the cytoplasm. Also there is a set of genes (7), related to ER that are also downregulated
#### on the total fraction.


########### Plot Total fractions RNA #######

mycolors <- c("blue", "red", "grey")
names(mycolors) <- c("DOWN", "UP", "NO")

#### Visulaizationn of downregulated genes in total fractions after 4 and 8 hours of auxin depletion ###############


p1<-ggplot(cur_LogMean_RNA_Data_total_4_select, aes(x = x, y = y,label=delabel,col=diffexpressed)) + 
  geom_point() +   xlab("Non treated") +
  ylab("4 hours Auxin treatment") + geom_text_repel() + scale_colour_manual(values = mycolors)

p2<-ggplot(cur_LogMean_RNA_Data_total_8_select, aes(x = x, y = y,label=delabel,col=diffexpressed)) + 
  geom_point() +   xlab("Non treated") +
  ylab("8 hours Auxin treatment") + geom_text_repel() + scale_colour_manual(values = mycolors)

#arrangment in a grid of the previous plots

layout <- rbind(c(1, 2))
grid.arrange(p1, p2,layout_matrix=layout)


#### We see how there is way less downregulation of genes in the total fraction compared to the ER fractions.
#### what can inform us about the fact that in general, depletion of SRP72 causes misloclization and a set of genes
#### related to the ER are downregulated significantly.




