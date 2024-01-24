library(DESeq2)
library(edgeR)
library(fdrtool)

library(riborex)
library(EnsDb.Hsapiens.v86)
library(dplyr)
library(biomaRt)
library(ggplot2)

total_fractions<- read.table(file="/Volumes/AG_Landthaler/amanukyan/Projects/RiboSeq/Pipeline/run/sofia_rna/output/count.txt", header = T, sep = "\t")
colnames(total_fractions) <- lapply(colnames(total_fractions), function(x) strsplit(x, split = "_S[0-9]+_")[[1]][1])
rownames(total_fractions) <- total_fractions[,1]
total_fractions <- total_fractions[,c(23:30)]
colnames(total_fractions) <- gsub("output.bam_files.", "", colnames(total_fractions))
colnames(total_fractions) <- gsub(".bam", "", colnames(total_fractions))


####
# DE analysis ####
####

# remove low count genes
count_gene <- apply(total_fractions, 1, max) 
rnadata <- total_fractions[which(count_gene > 10),]

rnaCond <- c("ctrl","ctrl","ctrl","time 4 hours","time 4 hours","time 4 hours","time 8 hours","time 8 hours")
rnaCond<-factor(rnaCond)
# get riborex data
print(head(rnadata))

# get regular deseq2 results

dds <- DESeqDataSetFromMatrix(countData = rnadata,
                              DataFrame(rnaCond),
                              design= ~ rnaCond)
dds <- DESeq(dds)
unique_rnaCond <- as.character(rev(unique(rnaCond)))
rna_results_total <- NULL
for(i in 1:(length(unique_rnaCond)-1)){
  for(j in (i+1):(length(unique_rnaCond))){
    comparison <- c("rnaCond", unique_rnaCond[i], unique_rnaCond[j])
    print(comparison)
    cur_results <- as.data.frame(results(dds, contrast = comparison))
    rna_results_total <- rbind(rna_results, data.frame(cur_results, gene = rownames(cur_results), comparison = paste0(comparison[2:3], collapse = "_vs_")))
  }
}
rna_results_total$padj[is.na(rna_results$padj)] <- 1
rna_results_total_8<-filter(rna_results_total, comparison == "time 8 hours_vs_ctrl")

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


datax <- data.frame(LogFC_RNA = rna_results_total_8$log2FoldChange, Sig_RNA = as.numeric(rna_results_total_8$padj), gene = rna_results_total_8$gene)
datax$Sig_RNA[is.na(datax$Sig_RNA)] <- 1
temp <- rep("NonSig", nrow(datax))
temp[datax$Sig_RNA < 0.05 & abs(datax$LogFC_RNA) > 1] <- "Sig"
datax$Sig <- factor(temp, levels = c("NonSig", "Sig"))
cur_LogMean_RNA_Data <- LogMean_RNA_Data[datax$gene,]

cur_LogMean_RNA_Data_total_8 <- data.frame(x = as.numeric(cur_LogMean_RNA_Data[,1]), 
                                     y = as.numeric(cur_LogMean_RNA_Data[,3]), 
                                     Sig = as.character(datax$Sig))

# add a column of NAs
cur_LogMean_RNA_Data_total_8$diffexpressed <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
cur_LogMean_RNA_Data_total_8$diffexpressed[datax$LogFC_RNA > 1.5] <- "UP"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
cur_LogMean_RNA_Data_total_8$diffexpressed[datax$LogFC_RNA < -1] <- "DOWN"

cur_LogMean_RNA_Data_total_8$gene_symbol<-datax$gene
cur_LogMean_RNA_Data_total_8$delabel <- NA
cur_LogMean_RNA_Data_total_8$delabel[cur_LogMean_RNA_Data_total_8$diffexpressed != "NO"] <- cur_LogMean_RNA_Data_total_8$gene_symbol[cur_LogMean_RNA_Data_total_8$diffexpressed != "NO"]

genes_down_RNA_total_8<-cur_LogMean_RNA_Data_total_8%>%
  dplyr::filter(diffexpressed == "DOWN")%>%
  dplyr::filter(Sig=="Sig")%>%
  pull(gene_symbol)



#list of genes with signal peptide (339 genes)

list_genesER<-read.table(file="/Volumes/AG_Landthaler/sofia/riboSRP72/sp_list.txt", header = T, sep = "\t")

genes_SP_down_total<-list_genesER[list_genesER$Symbol %in% genes_down_RNA_total_8,]

genes_SP_down_total<-print(genes_SP_down_total$Symbol)

sp_total<-enrichGO(genes_SP_down_total, OrgDb='org.Hs.eg.db', keyType= "SYMBOL",ont= "ALL",pAdjustMethod = "BH",pvalueCutoff  = 0.05,qvalueCutoff  = 0.05)

dotplot(sp_total)

#what we observe is that 297 from 339 genes with SP are downregulated in the total farction contrary to the ER fraction,
#when we apply sig in the filtering a bit up, the number goes down to 44.


