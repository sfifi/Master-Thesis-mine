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

####
# Import data ####
####

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

####
# Visualize ####
####

######RNA cond2 vs cond1, in this case I clarified above: 
#rna_results_8<-filter(rna_results, comparison == "time 8 hours_vs_ctrl") that I want to compare t8 vs ctrl

datax_4 <- data.frame(LogFC_RNA = rna_results_4$log2FoldChange, Sig_RNA = as.numeric(rna_results_4$padj), gene = rna_results_4$gene)
datax_4$Sig_RNA[is.na(datax_4$Sig_RNA)] <- 1
temp_4 <- rep("NonSig", nrow(datax_4))
temp_4[datax_4$Sig_RNA < 0.05 & abs(datax_4$LogFC_RNA) > 1] <- "Sig"
datax_4$Sig <- factor(temp_4, levels = c("NonSig", "Sig"))
cur_LogMean_RNA_Data_4 <- LogMean_RNA_Data[datax_4$gene,]


datax_8 <- data.frame(LogFC_RNA = rna_results_8$log2FoldChange, Sig_RNA = as.numeric(rna_results_8$padj), gene = rna_results_8$gene)
datax_8$Sig_RNA[is.na(datax_8$Sig_RNA)] <- 1
temp_8 <- rep("NonSig", nrow(datax_8))
temp_8[datax_8$Sig_RNA < 0.05 & abs(datax_8$LogFC_RNA) > 1] <- "Sig"
datax_8$Sig <- factor(temp_8, levels = c("NonSig", "Sig"))
cur_LogMean_RNA_Data_8 <- LogMean_RNA_Data[datax_8$gene,]


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
cur_LogMean_RNA_Data_4_select$diffexpressed[datax_4$LogFC_RNA < -1 & datax_8$Sig == "Sig"] <- "DOWN"
cur_LogMean_RNA_Data_4_select$gene_symbol<-datax_4$gene
cur_LogMean_RNA_Data_4_select$delabel <- NA
cur_LogMean_RNA_Data_4_select$delabel[cur_LogMean_RNA_Data_4_select$gene_symbol %in% individual_genes] <- cur_LogMean_RNA_Data_4_select$gene_symbol[cur_LogMean_RNA_Data_4_select$gene_symbol %in% individual_genes]

### For the condition of 8 hours

# add a column of NAs
cur_LogMean_RNA_Data_8_select$diffexpressed <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
cur_LogMean_RNA_Data_8_select$diffexpressed[datax_8$LogFC_RNA > 1 & datax_8$Sig == "Sig"] <- "UP"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
cur_LogMean_RNA_Data_8_select$diffexpressed[datax_8$LogFC_RNA < -1 & datax_8$Sig == "Sig"] <- "DOWN"
cur_LogMean_RNA_Data_8_select$gene_symbol<-datax_8$gene
#selection of UP genes with Sig

genes_up_RNA<-cur_LogMean_RNA_Data_8_select%>%
  dplyr::filter(diffexpressed == "UP")%>%
  dplyr::filter(Sig=="Sig")%>%
  pull(gene_symbol)
up_RNA_8_go<- enrichGO(genes_up_RNA, OrgDb='org.Hs.eg.db', keyType= "SYMBOL",ont= "ALL",pAdjustMethod = "BH",pvalueCutoff  = 0.05,qvalueCutoff  = 0.05)
up_RNA_8_go<-as.data.frame((up_RNA_8_go))
protein_folding_genes <- up_RNA_8_go$geneID[up_RNA_8_go$Description == "protein folding" |  up_RNA_8_go$Description == "response to unfolded protein" |
                                              up_RNA_8_go$Description == "response to topologically incorrect protein" |  up_RNA_8_go$Description == "'de novo' post-translational protein folding" |  up_RNA_8_go$Description == "response to heat"]
print(protein_folding_genes)
# Split the strings into a list of vectors
split_genes <- strsplit(protein_folding_genes, split = "/")
# Unlist the list to get a single vector of genes
individual_genes <- unlist(split_genes)
print(individual_genes)

#labelling
cur_LogMean_RNA_Data_8_select$delabel <- NA
cur_LogMean_RNA_Data_8_select$delabel[cur_LogMean_RNA_Data_8_select$gene_symbol %in% individual_genes] <- cur_LogMean_RNA_Data_8_select$gene_symbol[cur_LogMean_RNA_Data_8_select$gene_symbol %in% individual_genes]

mycolors <- c("blue", "red", "grey")
names(mycolors) <- c("DOWN", "UP", "NO")

#### Visulaizationn of downregulated genes after 4 and 8 hours of auxin depletion

p1 <- ggplot(cur_LogMean_RNA_Data_4_select, aes(x = x, y = y, label = delabel, col = diffexpressed)) + 
  geom_point() + 
  xlab("Non treated") +
  ylab("4 hours Auxin treatment") + 
  geom_text_repel() + 
  scale_colour_manual(values = mycolors) + 
  xlim(-5, 10) + 
  ylim(-6, 10) + 
  theme_minimal() +
  theme(
    axis.title.x = element_text(size = rel(2)), # make x axis label 1.5 times larger
    axis.title.y = element_text(size = rel(2)), # make y axis label 1.5 times larger
    plot.title = element_text(size = rel(2)),
    legend.text = element_text(size = rel(1)),
    legend.title = element_text(size = rel(1.5))# make title 2 times larger
  )  + geom_text_repel(max.overlaps=40,size=6,nudge_y = 2)

p2 <- ggplot(cur_LogMean_RNA_Data_8_select, aes(x = x, y = y, label = delabel, col = diffexpressed)) + 
  geom_point() + 
  xlab("Non treated") +
  ylab("8 hours Auxin treatment") + 
  geom_text_repel() + 
  scale_colour_manual(values = mycolors) + 
  xlim(-5, 10) + 
  ylim(-6, 10) + 
  theme_minimal() +
  theme(
    axis.title.x = element_text(size = rel(2)), # make x axis label 1.5 times larger
    axis.title.y = element_text(size = rel(2)), # make y axis label 1.5 times larger
    plot.title = element_text(size = rel(2)),
    legend.text = element_text(size = rel(1)),
    legend.title = element_text(size = rel(1.5))# make title 2 times larger
  ) + geom_text_repel(max.overlaps=40,size=6,nudge_y = 2)


#arrangment in a grid of the previous plots

layout <- rbind(c(1, 2))
grid.arrange(p1, p2,layout_matrix=layout)

####### working with RNA analyses #########

#which enrichment of genes are the ones up regulated?

dotplot(up_RNA_8_go)

RNA_up_8_gost<-gost(query = genes_up_RNA, 
                    organism = 'hsapiens',significant=T)
j<-gostplot(RNA_up_8_gost,capped = TRUE, interactive = FALSE)
IDs<-up_RNA_8_go@result[["ID"]]
select<-IDs[1:25]
publish_gostplot(j,highlight_terms = select)

######### what are the ones down regulated? ##########

genes_down_RNA_8<-cur_LogMean_RNA_Data_8_select%>%
  dplyr::filter(diffexpressed == "DOWN")%>%
  dplyr::filter(Sig=="Sig")%>%
  pull(gene_symbol)


genes_down_RNA_4<-cur_LogMean_RNA_Data_4_select%>%
  dplyr::filter(diffexpressed == "DOWN")%>%
  dplyr::filter(Sig=="Sig")%>%
  pull(gene_symbol)


#Visualization of enrichment of downregulated genes after 8 hours of SRP72 depletion.

down_RNA_8_go<- enrichGO(genes_down_RNA_8, OrgDb='org.Hs.eg.db', keyType= "SYMBOL",ont= "ALL",pAdjustMethod = "BH",pvalueCutoff  = 0.05,qvalueCutoff  = 0.05)

dotplot(down_RNA_8_go)

RNA_down_8_gost<-gost(query = genes_down_RNA, 
                      organism = 'hsapiens',significant=FALSE)

j<-gostplot(RNA_down_8_gost,capped = TRUE, interactive = FALSE)

IDs<-down_RNA_8_go@result[["ID"]]

select_down<-IDs[1:10]

publish_gostplot(j,highlight_terms = select_down)

down_RNA_4_go<- enrichGO(genes_down_RNA_4, OrgDb='org.Hs.eg.db', keyType= "SYMBOL",ont= "ALL",pAdjustMethod = "BH",pvalueCutoff  = 0.05,qvalueCutoff  = 0.05)

dotplot(down_RNA_4_go)



# TE vs GEX LFC

common_genes <- intersect(rna_results_8$gene, ribo_results_8$gene)
rnavis <- rna_results_8[match(common_genes, rna_results_8$gene),]
ribovis <- ribo_results_8[match(common_genes, ribo_results_8$gene),]
datax_TE <- data.frame(logFC_TE = ribovis$log2FoldChange, LogFC_RNA = rnavis$log2FoldChange, Sig_Ribo = ribovis$padj, Sig_RNA = rnavis$padj)
datax_TE$Sig_Ribo[is.na(datax$Sig_Ribo)] <- 1
temp <- rep("NonSig", nrow(datax_TE))
temp[datax_TE$Sig_Ribo < 0.05 & abs(datax_TE$logFC_TE) > 1.5] <- "Sig"
datax_TE$Sig_TE <- factor(temp, levels = c("NonSig", "Sig"))

temp <- rep("NonSig", nrow(datax_TE))
temp[datax$Sig_RNA < 0.05 & abs(datax_TE$LogFC_RNA) > 1.5] <- "Sig"
datax_TE$Sig_GEX <- factor(temp, levels = c("NonSig", "Sig"))

Sig_all <- rep("NonSig", nrow(datax_TE))
Sig_all[datax_TE$Sig_TE == "Sig" & datax_TE$Sig_GEX == "Sig"] <- "Sig (TE/GEX)"
Sig_all[datax_TE$Sig_TE == "Sig" & datax_TE$Sig_GEX == "NonSig"] <- "Sig (TE)"
Sig_all[datax_TE$Sig_TE == "NonSig" & datax_TE$Sig_GEX == "Sig"] <- "Sig (GEX)"
datax_TE$Sig_all <- factor(Sig_all, levels = rev(c("NonSig", "Sig (GEX)", "Sig (TE)", "Sig (TE/GEX)")))

# add a column of NAs
datax_TE$diffexpressed <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
datax_TE$diffexpressed[datax_TE$LogFC_RNA > 1.5] <- "UP"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
datax_TE$diffexpressed[datax_TE$LogFC_RNA < -2] <- "DOWN"

datax_TE$gene_symbol<-common_genes
datax_TE$delabel <- NA
datax_TE$delabel[datax_TE$Sig_all != "NonSig"] <- datax_TE$gene_symbol[datax_TE$Sig_all != "NonSig"]


mycolors <- c("grey","blue", "red", "darkgreen")
names(mycolors) <- c("NonSig","Sig (TE/GEX)", "Sig (TE)", "Sig (GEX)")


ggplot(datax_TE, aes(x = LogFC_RNA, y = logFC_TE, label=delabel,col=Sig_all)) + 
  geom_point() + labs(title ="TE vs GEX") + geom_text_repel() + scale_colour_manual(values = mycolors)


#enrichment of pathways for TE/GEX

select_sig_TEGEX<-datax_TE[datax_TE$Sig_all == "Sig (TE/GEX)",]

genes_TEGEX<-select_sig_TEGEX[,9]

go<- enrichGO(genes_TEGEX, OrgDb='org.Hs.eg.db', keyType= "SYMBOL",ont= "ALL",pAdjustMethod = "BH",pvalueCutoff  = 0.05,qvalueCutoff  = 0.05)

#for gprofiler2:

TEGEXgost<-gost(query = genes_TEGEX, 
                organism = 'hsapiens',significant=FALSE)

gostplot(TEGEXgost,capped = TRUE, interactive = TRUE)


#for clusterprofiler:

entrez_ids <- mapIds(org.Hs.eg.db, keys = genes_TEGEX,column="ENTREZID", keytype = "SYMBOL",select=genes_TEGEX)

TEGEXgsego<-gseGO(entrez_ids, OrgDb="org.Hs.eg.db",ont= "ALL",decreasing=TRUE,eyType = "ENTREZID",exponent = 1,minGSSize = 3,maxGSSize = 500,)

dotplot(TEGEXgsego,showCategory=30)

#enrichment of pathways for GEX

e<-datax_TE[datax_TE$Sig_GEX == "Sig",]
f<-e[e$diffexpressed=="DOWN",]
g<-f[,9]


go2<- enrichGO(g, OrgDb='org.Hs.eg.db', keyType= "SYMBOL",ont= "ALL",pAdjustMethod = "BH",pvalueCutoff  = 0.05,qvalueCutoff  = 0.05)

GEXgost<-gost(query =g, organism = "hsapiens", significant=FALSE)

gostplot(GEXgost,capped = TRUE, interactive = TRUE)


##enrichment of pathways for TE

Sig_TE<-datax_TE[datax_TE$Sig_all == "Sig (TE)",]

#Significant genes TE Logfold down

Sig_TE_down<-Sig_TE[Sig_TE$logFC_TE < -1.5,]
Sig_TE_genes_down<-Sig_TE_down[,9]

##Significant genes TE Logfold up

Sig_TE_up<-Sig_TE[Sig_TE$logFC_TE > 1.5,]
Sig_TE_genes_up<-Sig_TE_up[,9]


multi_TE<-gost(list("up-regulated"=Sig_TE_genes_up,"down-regulated"=Sig_TE_genes_down), organism = "hsapiens", multi_query=F,evcodes=T, significant=T)

multi<-gostplot(multi_TE,capped = TRUE, interactive = F)

genes_TE_up<-enrichGO(Sig_TE_genes_up, OrgDb='org.Hs.eg.db', keyType= "SYMBOL",ont= "ALL")

dotplot(genes_TE_up,showCategory=30)

IDs<-genes_TE_up@result[["ID"]]

publish_gostplot(multi,highlight_terms = IDs)


# modify the g:Profiler data frame
# modify the g:Profiler data frame
gp_mod <- multi_TE$result[,c("query", "source", "term_id",
                             "term_name", "p_value", "query_size",
                             "intersection_size", "term_size",
                             "effective_domain_size", "intersection")]
gp_mod$GeneRatio <- paste0(gp_mod$intersection_size,  "/", gp_mod$query_size)
gp_mod$BgRatio <- paste0(gp_mod$term_size, "/", gp_mod$effective_domain_size)
names(gp_mod) <- c("Cluster", "Category", "ID", "Description", "p.adjust",
                   "query_size", "Count", "term_size", "effective_domain_size",
                   "geneID", "GeneRatio", "BgRatio")
gp_mod$geneID <- gsub(",", "/", gp_mod$geneID)
row.names(gp_mod) <- gp_mod$ID
# define as compareClusterResult object
gp_mod_cluster <- new("compareClusterResult", compareClusterResult = gp_mod)
# define as enrichResult object
gp_mod_enrich <- new("enrichResult", result = gp_mod)

dotplot(gp_mod_cluster)

barplot(gp_mod_enrich, showCategory = 60, font.size = 16) +
  ggplot2::facet_grid(~Cluster) +
  ggplot2::ylab("Intersection size")


######## TE cond2 vs cond1 ? (exercise)


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

list_genes<-c("UPF1","ATF4","ATF5","PRPF6","PRPF19","SNRPN200","EIF3B","EIF1AD","NOP9","NOP14","WFS1","RNF40","USP19","HSF1")

cur_LogMean_TE_Data_8_select$delabel[cur_LogMean_TE_Data_8_select$gene_symbol %in% list_genes] <- cur_LogMean_TE_Data_8_select$gene_symbol[cur_LogMean_TE_Data_8_select$gene_symbol %in% list_genes]


mycolors <- c("blue", "red", "grey")
names(mycolors) <- c("DOWN", "UP", "NO")


p3 <- ggplot2::ggplot(cur_LogMean_TE_Data_4_select, aes(x = x, y = y, label = delabel, col = diffexpressed)) + 
  geom_point() + 
  xlab("Non treated") +
  ylab("4 hours Auxin treatment") + 
  geom_text_repel(max.overlaps = 40, size = 3, nudge_y = 2) + 
  scale_colour_manual(values = mycolors) + 
  xlim(-6, 10) + 
  ylim(-6, 10) + 
  theme_minimal() + 
  theme(
    axis.title.x = element_text(size = rel(2)), # make x axis label 2 times larger
    axis.title.y = element_text(size = rel(2)), # make y axis label 2 times larger
    plot.title = element_text(size = rel(2)),   # make title 2 times larger
    legend.text = element_text(size = rel(1)),  # make legend text 1 time larger
    legend.title = element_text(size = rel(1.5))# make legend title 1.5 times larger
  )

p4<-ggplot2::ggplot(cur_LogMean_TE_Data_8_select, aes(x = x, y = y,label=delabel,col=diffexpressed)) + 
  geom_point() +   xlab("Non treated") +
  ylab("8 hours Auxin treatment") + geom_text_repel(max.overlaps=40,size=6,nudge_y = 2) + scale_colour_manual(values = mycolors) + xlim(-6,10) + ylim(-6,10) + theme_minimal() +
  theme(
    axis.title.x = element_text(size = rel(2)), # make x axis label 1.5 times larger
    axis.title.y = element_text(size = rel(2)), # make y axis label 1.5 times larger
    plot.title = element_text(size = rel(2)),
    legend.text = element_text(size = rel(1)),
    legend.title = element_text(size = rel(1.5))# make title 2 times larger
  )


layout <- rbind(c(1, 2))
grid.arrange(p3, p4,layout_matrix=layout)

########go enrichment and gostplot#########

####for condition 4 hours

TE_4_up<-cur_LogMean_TE_Data_4_select[cur_LogMean_TE_Data_4_select$diffexpressed == "UP",]
TE_4_up<-TE_4_up[TE_4_up$Sig =="Sig",]
TE_4_up_genes<-TE_4_up[,5]

TE_4<-gost(TE_4_up_genes, organism = "hsapiens",significant=FALSE)

multi_4_up <- gostplot(TE_4,capped = TRUE, interactive = F)

genes_TE_4up<-enrichGO(TE_4_up_genes, OrgDb='org.Hs.eg.db', keyType= "SYMBOL",ont= "ALL")

IDs4<-genes_TE_4up@result[["ID"]]

publish_gostplot(multi_4_up,highlight_terms = IDs)

#####for condition 8 hours


TE_8_up<-cur_LogMean_TE_Data_8_select[cur_LogMean_TE_Data_8_select$diffexpressed == "UP",]
TE_8_up<-TE_8_up[TE_8_up$Sig=="Sig",]
TE_8_up_genes<-TE_8_up[,5]

TE_8<-gost(TE_8_up_genes, organism = "hsapiens",significant=T)

multi_8_up<-gostplot(TE_8,capped = TRUE, interactive = T)

genes_TE_8up<-enrichGO(TE_8_up_genes, OrgDb='org.Hs.eg.db', keyType= "SYMBOL",ont= "MF",pAdjustMethod = "BH",pvalueCutoff  = 0.05,qvalueCutoff  = 0.05)

dotplot(genes_TE_8up)
barplot(genes_TE_8up)

IDs8<-genes_TE_8up@result[["ID"]]

publish_gostplot(multi_8_up,highlight_terms = IDs)

##### TE down ######

TE_8_down<-cur_LogMean_TE_Data_8_select[cur_LogMean_TE_Data_8_select$diffexpressed == "DOWN",]
TE_8_down<-TE_8_down[TE_8_down$Sig=="Sig",]
TE_8_down_genes<-TE_8_down[,5]

TE_8_down<-gost(TE_8_down_genes, organism = "hsapiens",significant=T)

multi_8_down<-gostplot(TE_8_down,capped = TRUE, interactive = T)


genes_TE_8down<-enrichGO(TE_8_down_genes, OrgDb='org.Hs.eg.db', keyType= "SYMBOL",ont= "ALL",pAdjustMethod = "BH",pvalueCutoff  = 0.05,qvalueCutoff  = 0.05)

barplot(genes_TE_8down)


datax_ribo_8$diffexpressed <- "NO"

datax_ribo_8$diffexpressed[datax_ribo_8$LogFC_RNA > 1 & datax_ribo_8$Sig == "Sig"] <- "UP"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
datax_ribo_8$diffexpressed[datax_ribo_8$LogFC_RNA < -1 & datax_ribo_8$Sig == "Sig"] <- "DOWN"


d<-datax_ribo_8[datax_ribo_8$diffexpressed == "DOWN" & datax_ribo_8$Sig=="Sig",]


d2<-d%>%
  arrange(LogFC_RNA)

####### SP and TE ###########

list_genesER<-read.table(file="/Volumes/AG_Landthaler/sofia/riboSRP72/sp_list.txt", header = T, sep = "\t")

genes_SP_TE_8<-list_genesER[list_genesER$Symbol %in% TE_8_up$gene_symbol,]

genes_SP_TE_8<-print(genes_SP_TE_8$Symbol)


# genes PTPRF and AGPAT2 are significantly more translated after 8 h containing SP signal.


#Barplots genes

rnadata_tpm$gene<-rownames(rnadata_tpm)
colnames(rnadata_tpm)<-c("E0R1","E0R2","E4R1","E4R2","E8R1","E8R2", "gene")
ribodata_tpm$gene<-rownames(ribodata_tpm)
colnames(ribodata_tpm)<-c("E0R1","E0R2","E4R1","E4R2","E8R1","E8R2", "gene")


melted_rnatpm<-melt(rnadata_tpm, id.vars=c("gene"))
melted_rnatpm$condition <- gsub("R[0-9]+", "", melted_rnatpm$variable)
melted_ribotpm<-melt(ribodata_tpm, id.vars=c("gene"))
melted_ribotpm$condition <- gsub("R[0-9]+", "", melted_ribotpm$variable)

tablerna_tpm<-melted_rnatpm %>%
  group_by(condition,gene)%>%
  summarize(mean_value = mean(value))
tablerna_tpm$mean_value <- log(tablerna_tpm$mean_value + 1)

tableribo_tpm<-melted_ribotpm %>%
  group_by(condition,gene)%>%
  summarize(mean_value = mean(value))
tableribo_tpm$mean_value <- log(tableribo_tpm$mean_value + 1)

table_tpm_bar<-merge(tablerna_tpm,tableribo_tpm, by = c("gene", "condition"), all = TRUE)

#value x is from rna and y from ribo

grouped_table<-melt(table_tpm_bar,id.vars=c("gene","condition"))

desired_levels <- c("mean_value.x", "mean_value.y")

# grouped_table$variable <- factor(grouped_table$variable, levels = desired_levels)

grouped_table <- grouped_table %>%
  mutate(variable = recode(variable, 
                           "mean_value.x" = "rna", 
                           "mean_value.y" = "ribo"))

gene_put<-grouped_table[grouped_table$gene %in% c("GAPDH","EPPK1","H2BC11","H4C2","COL12A1","TIMM10"),]

ggplot(data = gene_put, aes(x = condition, y = value,fill=variable)) +
  geom_bar(stat = "identity",position="dodge") + labs(title = "RNA vs Ribo") +facet_wrap(~ gene)

##########corrplot ##########

library(debrowser)

normdata_ribo_ER <- getNormalizedMatrix(ribodata_tpm, method = "TMM")

library(corrplot)
plot<-all2all(normdata_ribo_ER,cex=2)

ribodata_tpm$gene<-rownames(ribodata_tpm)
ribodata$gene<-rownames(ribodata)
