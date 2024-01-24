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
# loading the additional packages
library(clusterProfiler)
library(enrichplot)
library(DOSE) # needed to convert to enrichResult object
library(gridExtra)
library(ggrepel)
library(reshape2)
library(tidyverse)


# library(riborex)
source("/Volumes/AG_Landthaler/sofia/DESeq2Rex.R")

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
rnadata <- rnadata[,c(22:23,25:26,28:29)]
colnames(rnadata) <- gsub("output.bam_files.", "", colnames(rnadata))


ribodata <- read.table(file="/Volumes/AG_Landthaler/amanukyan/Projects/RiboSeq/Pipeline/run/sofia/output/FeatureCounts/cds_counts.txt", header = T, sep = "\t")
colnames(ribodata) <- lapply(colnames(ribodata), function(x) strsplit(x, split = "_S[0-9]+_")[[1]][1])
rownames(ribodata) <- ribodata[,1]
ribodata <- ribodata[,-1]
ribodata <- ribodata[,c(20:21,23:26)]
colnames(ribodata) <- gsub("output.bam_files.", "", colnames(ribodata))
colnames(ribodata) <- gsub(".bam", "", colnames(ribodata))

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
# Calculate TE for CDS and UTR ####
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

# calculate log TE mean for each condition. Here is calculated the values of TE.

cond_names <- riboCond
common_genes <- intersect(rownames(ribodata_tpm), rownames(rnadata_tpm))

#calculation of TE based on TPM from RNA and ribo

TE_Data <- ribodata_tpm[common_genes,]/rnadata_tpm[common_genes,]
Mean_TE_Data <- aggregate(t(TE_Data), list(cond_names), mean)
Mean_TE_Data <- t(Mean_TE_Data)
colnames(Mean_TE_Data) <- Mean_TE_Data[1,]
Mean_TE_Data <- data.frame(apply(Mean_TE_Data[-1,],2,as.numeric), row.names = rownames(Mean_TE_Data[-1,]))
LogMean_TE_Data <- log(Mean_TE_Data + 0.001)

#combining TPM with significancy from Deseq for both conditions 4 and 8 hours


datax_ribo_4 <- data.frame(LogFC_RNA = ribo_results_4$log2FoldChange, Sig_RNA = as.numeric(ribo_results_4$padj), gene = ribo_results_4$gene)
datax_ribo_4$Sig_RNA[is.na(datax_ribo_4$Sig_RNA)] <- 1
temp_4 <- rep("NonSig", nrow(datax_ribo_4))
temp_4[datax_ribo_4$Sig_RNA < 0.05 & abs(datax_ribo_4$LogFC_RNA) > 1] <- "Sig"
datax_ribo_4$Sig <- factor(temp_4, levels = c("NonSig", "Sig"))
cur_LogMean_TE_Data_4 <- LogMean_TE_Data[datax_ribo_4$gene,]

datax_ribo_8 <- data.frame(LogFC_RNA = ribo_results_8$log2FoldChange, Sig_RNA = as.numeric(ribo_results_8$padj), gene = ribo_results_8$gene)
datax_ribo_8$Sig_RNA[is.na(datax_ribo_8$Sig_RNA)] <- 1
temp_8 <- rep("NonSig", nrow(datax_ribo_8))
temp_8[datax_ribo_8$Sig_RNA < 0.05 & abs(datax_ribo_8$LogFC_RNA) > 1] <- "Sig"
datax_ribo_8$Sig <- factor(temp_8, levels = c("NonSig", "Sig"))
cur_LogMean_TE_Data_8 <- LogMean_TE_Data[datax_ribo_8$gene,]

#x and y are normalized data, x=ctrl and y= 4 h or 8 (column 2 or 3 respectively)

cur_LogMean_TE_Data_4_select <- data.frame(x = as.numeric(cur_LogMean_TE_Data_4[,1]), 
                                           y = as.numeric(cur_LogMean_TE_Data_4[,2]), 
                                           Sig = as.character(datax_ribo_4$Sig))


cur_LogMean_TE_Data_8_select <- data.frame(x = as.numeric(cur_LogMean_TE_Data_8[,1]), 
                                           y = as.numeric(cur_LogMean_TE_Data_8[,3]), 
                                           Sig = as.character(datax_ribo_8$Sig))

# add a column of NAs
cur_LogMean_TE_Data_4_select$diffexpressed <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
cur_LogMean_TE_Data_4_select$diffexpressed[datax_ribo_4$LogFC_RNA > 1 & datax_ribo_4$Sig == "Sig"] <- "UP"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
cur_LogMean_TE_Data_4_select$diffexpressed[datax_ribo_4$LogFC_RNA < -1 & datax_ribo_4$Sig == "Sig"] <- "DOWN"

cur_LogMean_TE_Data_4_select$gene_symbol<-datax_ribo_4$gene
cur_LogMean_TE_Data_4_select$delabel <- NA
cur_LogMean_TE_Data_4_select$delabel[cur_LogMean_TE_Data_4_select$diffexpressed != "NO"] <- cur_LogMean_TE_Data_4_select$gene_symbol[cur_LogMean_TE_Data_4_select$diffexpressed != "NO"]

# add a column of NAs
cur_LogMean_TE_Data_8_select$diffexpressed <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
cur_LogMean_TE_Data_8_select$diffexpressed[datax_ribo_8$LogFC_RNA > 1 & datax_ribo_8$Sig == "Sig"] <- "UP"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
cur_LogMean_TE_Data_8_select$diffexpressed[datax_ribo_8$LogFC_RNA < -1 & datax_ribo_8$Sig == "Sig"] <- "DOWN"

cur_LogMean_TE_Data_8_select$gene_symbol<-datax_ribo_8$gene
cur_LogMean_TE_Data_8_select$delabel <- NA
cur_LogMean_TE_Data_8_select$delabel[cur_LogMean_TE_Data_8_select$diffexpressed != "NO"] <- cur_LogMean_TE_Data_8_select$gene_symbol[cur_LogMean_TE_Data_8_select$diffexpressed != "NO"]


mycolors <- c("blue", "red", "grey")
names(mycolors) <- c("DOWN", "UP", "NO")


p3<-ggplot(cur_LogMean_TE_Data_4_select, aes(x = x, y = y,label=delabel,col=diffexpressed)) + xlim(-6,10) + ylim(-6,10) + theme_minimal() +
  geom_point() +   xlab("Non treated") +
  ylab("4 hours Auxin treatment") + geom_text_repel() + scale_colour_manual(values = mycolors)+
  theme(
    axis.title.x = element_text(size = rel(2)), # make x axis label 1.5 times larger
    axis.title.y = element_text(size = rel(2)), # make y axis label 1.5 times larger
    plot.title = element_text(size = rel(2)),
    legend.text = element_text(size = rel(1)),
    legend.title = element_text(size = rel(1.5))# make title 2 times larger
  )

p4<-ggplot(cur_LogMean_TE_Data_8_select, aes(x = x, y = y,label=delabel,col=diffexpressed)) + 
  geom_point() +   xlab("Non treated") +
  ylab("8 hours Auxin treatment") + geom_text_repel() + scale_colour_manual(values = mycolors) + xlim(-6,10) + ylim(-6,10) + theme_minimal()+
  theme(
    axis.title.x = element_text(size = rel(2)), # make x axis label 1.5 times larger
    axis.title.y = element_text(size = rel(2)), # make y axis label 1.5 times larger
    plot.title = element_text(size = rel(2)),
    legend.text = element_text(size = rel(1)),
    legend.title = element_text(size = rel(1.5))# make title 2 times larger
  )


layout <- rbind(c(1, 2))
grid.arrange(p3, p4,layout_matrix=layout)

up_TE_total<-datax_ribo_8[datax_ribo_8$Sig=="Sig",]
sorted_genes_asc <- up_TE_total$gene[order(desc(up_TE_total$LogFC_RNA))]
sorted_genes_asc<-sorted_genes_asc[3:50]

TE_total_enrich_up<-enrichGO(sorted_genes_asc, OrgDb='org.Hs.eg.db', keyType= "SYMBOL",ont= "ALL",pAdjustMethod = "BH",pvalueCutoff  = 0.05,qvalueCutoff  = 0.05)
dotplot(TE_total_enrich_up)

sorted_genes_desc <- up_TE_total$gene[order(up_TE_total$LogFC_RNA)]
sorted_genes_desc<-sorted_genes_desc[1:20]

TE_total_enrich_down<-enrichGO(sorted_genes_desc, OrgDb='org.Hs.eg.db', keyType= "SYMBOL",ont= "ALL",pAdjustMethod = "BH",pvalueCutoff  = 0.05,qvalueCutoff  = 0.05)
dotplot(TE_total_enrich_down)
