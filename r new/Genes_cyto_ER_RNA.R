library(DESeq2)
library(EnsDb.Hsapiens.v86)
library(dplyr)
library(biomaRt)
library(ggplot2)
library(readxl)

setwd("/Users/sgarcia/R")

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

library(readxl)
library(tidyverse)

C_E <- read_excel("/Users/sgarcia/Downloads/cyto_ER_genes.xls",sheet = 1,col_names = T, col_types = "text")
C_genes <- C_E$Symbol[C_E$localization=="cytosolic"]
E_genes <- C_E$Symbol[C_E$localization=="membrane"]

#datax_cyto_genes <- rna_results[rna_results$log2FoldChange < -0.5,]
#datax <- rna_results[rna_results$log2FoldChange > 1,]

######## HISTOGRAM PLOT ER CYTO ############


ggplot2::ggplot(datax_CE, aes(log2FoldChange, fill=localisation))+geom_histogram(bins=250)+
  scale_fill_manual(values=c("dodgerblue4","orange2","grey")) + theme_minimal()+
  theme(axis.text = element_text(size = 16))

cytosolic_genes <- datax_CE[datax_CE$localisation == "cytosolic",]
#cytosolic_genes<-write.table(df, file = "cytosolic_genes.txt", sep = "\t", row.names = FALSE)
cytosolic_genes <- cytosolic_genes$gene

membrane_genes <- datax_CE[datax_CE$localisation=="membrane",]
#membrane_genes<-write.table(df, file = "membrane_genes.txt", sep = "\t", row.names = FALSE)
membrane_genes <- membrane_genes$gene

undefined_genes <- datax_CE[datax_CE$localisation=="undefined",]
undefined_genes <- undefined_genes$gene

###### get the byotype of these genes #######

ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
genes <- datax_CE$gene 

gene_info <- getBM(
  attributes = c("external_gene_name", "gene_biotype"),
  filters = "external_gene_name",
  values = genes,
  mart = ensembl
)

datax_CE <- datax_CE %>%
  left_join(gene_info, by = c("gene" = "external_gene_name"))%>%
  mutate(gene_biotype = coalesce(gene_biotype, "noinfo"))

membrane<-datax_CE[datax_CE$localisation=="membrane",]

# Count the occurrences of each gene biotype
gene_counts <- table(membrane$gene_biotype)

# Print the counts
print(gene_counts)

############################## fromSPanaylses script #####################

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

datax_8 <- data.frame(LogFC_RNA = rna_results_8$log2FoldChange, Sig_RNA = as.numeric(rna_results_8$padj), gene = rna_results_8$gene)
datax_8$Sig_RNA[is.na(datax_8$Sig_RNA)] <- 1
temp_8 <- rep("NonSig", nrow(datax_8))
temp_8[datax_8$Sig_RNA < 0.05 & abs(datax_8$LogFC_RNA) > 1] <- "Sig"
datax_8$Sig <- factor(temp_8, levels = c("NonSig", "Sig"))

# add a column of NAs
datax_8$diffexpressed <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
datax_8$diffexpressed[datax_8$LogFC_RNA > 1 & datax_8$Sig == "Sig"] <- "UP"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
datax_8$diffexpressed[datax_8$LogFC_RNA < -1 & datax_8$Sig == "Sig"] <- "DOWN"


############# SET separation ####################

genes_MITO<-datax_8$gene %in% result_df$Gene_Symbol
genes_SP<-datax_8$gene %in% list_genesSP$Symbol
genes_TM<-datax_8$gene %in% list_genesTM$Symbol

datax_8$Set <- NA
datax_8$Set[genes_MITO] <- "mtRNA"
datax_8$Set[genes_SP] <- "SP"
datax_8$Set[genes_TM] <- "TM"
categories_to_exclude <- c("mtRNA", "SP", "TM")
datax_8$Set <- ifelse(datax_8$Set %in% categories_to_exclude, datax_8$Set, "No SP/TM or mtRNA")

datax_8<-datax_8[datax_8$Sig=="Sig",]
datax_8$localisation <- NA
datax_8$localisation[datax_8$gene %in% membrane_genes] <- "membrane"
datax_8$localisation[datax_8$gene %in% cytosolic_genes] <- "cytosolic"
datax_8$localisation[datax_8$gene %in% undefined_genes] <- "undefined"

datax_4 <- data.frame(LogFC_RNA = rna_results_4$log2FoldChange, Sig_RNA = as.numeric(rna_results_4$padj), gene = rna_results_4$gene)
datax_4$Sig_RNA[is.na(datax_4$Sig_RNA)] <- 1
temp_4 <- rep("NonSig", nrow(datax_4))
temp_4[datax_4$Sig_RNA < 0.05 & abs(datax_4$LogFC_RNA) > 1] <- "Sig"
datax_4$Sig <- factor(temp_4, levels = c("NonSig", "Sig"))

# add a column of NAs
datax_4$diffexpressed <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
datax_4$diffexpressed[datax_4$LogFC_RNA > 1 & datax_4$Sig == "Sig"] <- "UP"
####################
##################
##################
#############
##################
################
###############

datax_8_up <- datax_8[datax_8$diffexpressed=="UP",]
datax_8_down <- datax_8[datax_8$diffexpressed=="DOWN",]
common_genes_up <- datax_8_up[datax_8_up$gene %in% datax_4$gene,]
common_genes_down <- datax_8_down[datax_8_down$gene %in% datax_4$gene,]
datax_4$diffexpressed[datax_4$gene %in% common_genes_up$gene] <- "UP"
datax_4$diffexpressed[datax_4$gene %in% common_genes_down$gene] <- "DOWN"
common_genes <- datax_8[datax_8$gene %in% datax_4$gene,]
datax_4 <- datax_4[datax_4$gene %in% common_genes$gene,]

############# SET separation ####################

genes_MITO<-datax_4$gene %in% result_df$Gene_Symbol
genes_SP<-datax_4$gene %in% list_genesSP$Symbol
genes_TM<-datax_4$gene %in% list_genesTM$Symbol
datax_4$Set <- NA
datax_4$Set[genes_MITO] <- "mtRNA"
datax_4$Set[genes_SP] <- "SP"
datax_4$Set[genes_TM] <- "TM"
categories_to_exclude <- c("mtRNA", "SP", "TM")
datax_4$Set <- ifelse(datax_4$Set %in% categories_to_exclude, datax_4$Set, "No SP/TM or mtRNA")

datax_4$localisation <- NA
datax_4$localisation[datax_4$gene %in% membrane_genes] <- "membrane"
datax_4$localisation[datax_4$gene %in% cytosolic_genes] <- "cytosolic"
datax_4$localisation[datax_4$gene %in% undefined_genes] <- "undefined"

###### addition of byotype ###########

genes_4 <- datax_4$gene 
gene_info_4 <- getBM(
  attributes = c("external_gene_name", "gene_biotype"),
  filters = "external_gene_name",
  values = genes_4,
  mart = ensembl
)
datax_4 <- datax_4 %>%
  left_join(gene_info_4, by = c("gene" = "external_gene_name"))%>%
  mutate(gene_biotype = coalesce(gene_biotype, "noinfo"))

genes_8 <- datax_8$gene 
gene_info_8 <- getBM(
  attributes = c("external_gene_name", "gene_biotype"),
  filters = "external_gene_name",
  values = genes_8,
  mart = ensembl
)
datax_8 <- datax_8 %>%
  left_join(gene_info_8, by = c("gene" = "external_gene_name"))%>%
  mutate(gene_biotype = coalesce(gene_biotype, "noinfo"))

######## Source addition

# Add a column to datax_8 to indicate the source
datax_8 <- datax_8 %>%
  mutate(Source = "datax_8")
# Add a column to datax_4 to indicate the source
datax_4 <- datax_4 %>%
  mutate(Source = "datax_4")

# Combine the two data frames into one
combined_data <- rbind(datax_8, datax_4)
combined_data<- combined_data %>%
  mutate(Source = ifelse(Source == "datax_4", "4 h", ifelse(Source == "datax_8", "8 h", Source)))
combined_data_plot <-  combined_data
combined_data_plot$Group <- paste0(combined_data_plot$Set, ' (', combined_data_plot$Source, ')')
combined_data_plot <- combined_data_plot %>%
  filter(localisation %in% c("cytosolic", "membrane") & localisation != "undefined")
count_data <- combined_data_plot %>%
  group_by(Group,localisation) %>%
  summarise(Count = n(), .groups = 'drop')

# Assuming melted_data is your data frame
# Calculate the count of observations for each group

# Create the violin plot UP

ggplot(combined_data_plot, aes(x = Group, y = LogFC_RNA)) +
  geom_violin(aes(fill = Set), alpha = 0.5) + 
  geom_line(aes(group = gene), linewidth = 0.4, color = "gray") +
  geom_hline(yintercept = 0) + 
  geom_jitter(data = combined_data_plot, aes(x = Group, y = LogFC_RNA), width = .1, alpha = .3) +
  labs(title = "LogFC RNA on ER", x = "Sets", y = "LogFC RNA") +
  scale_fill_manual(values = c("SP" = "orange", "TM" = "dark green", "mtRNA" = "purple", "No SP/TM or mtRNA" = "grey")) +
  theme_minimal() +
  geom_text(data = count_data, aes(x = Group, y = Inf, label = Count), vjust = 1.5, size = 3, position = position_jitter(width = 0.1)) + 
  ylim(-5, 5) + 
  facet_grid(~localisation) +
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

####### simplified plot ##########

# Filter the data to include only "cytosolic" and "membrane" localizations
combined_data <- rbind(datax_8, datax_4)
combined_data<- combined_data %>%
  mutate(Source = ifelse(Source == "datax_4", "4 hours", ifelse(Source == "datax_8", "8 hours", Source)))
combined_data_plot <-  combined_data
filtered_data <- combined_data_plot %>%
  filter(localisation %in% c("cytosolic", "membrane") & localisation != "undefined")
count_data <- filtered_data %>%
  group_by(localisation,Source) %>%
  summarise(Count = n(), .groups = 'drop')


# Create the plot with the filtered data without using "Group"
ggplot(filtered_data, aes(x = Source, y = LogFC_RNA)) +
  geom_violin(aes(fill = localisation), alpha = 0.5) + 
  geom_line(aes(group = gene), linewidth = 0.4, color = "gray") +
  geom_hline(yintercept = 0) + 
  geom_jitter(data = filtered_data, aes(x = Source, y = LogFC_RNA), width = .1, alpha = .3) +
  labs(title = "LogFC RNA  on ER", x = "Time points", y = "LogFC RNA") +
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














######### enrich plots

###### Question: how many genes are cytosolic and membrane in my data and papers' data? ######
membrane_genes_down_mine <- datax_8_down[datax_8_down$gene %in% membrane_genes,]
cytosolic_genes_down_mine <- datax_8_down[datax_8_down$gene %in% cytosolic_genes,]
undefined_genes_down_mine <-datax_8_down[datax_8_down$gene %in% undefined_genes,]
C_genes_paper <- datax_8_down[datax_8_down$gene %in% C_genes,]
M_genes_paper <- datax_8_down[datax_8_down$gene %in% E_genes,]
nrow(membrane_genes_down_mine)
nrow(cytosolic_genes_down_mine)
nrow(C_genes_paper)
nrow(M_genes_paper)

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

######### There is a downregulation of genes supposed to be cytosolic ########
####### next, how many of those genes have SP or TM? #######
############# 37 genes from those 212 contain SP signal #######

SP_membrane_down <- list_genesSP[list_genesSP$Symbol %in% membrane_genes_down_mine$gene,]
nrow(SP_membrane_down)

TM_membrane_down <- list_genesTM[list_genesTM$Symbol %in% membrane_genes_down_mine$gene,]
nrow(TM_membrane_down)

MITO_membrane_down <- result_df[result_df$Gene_Symbol %in% membrane_genes_down_mine$gene,]
nrow(MITO_membrane_down)

######### Only 28 from those 212 contain TM domain #######
# what are those 100 genes then?


hundred_genes <- membrane_genes_down_mine[!(membrane_genes_down_mine$gene %in% SP_membrane_down) & !(membrane_genes_down_mine$gene %in% TM_membrane_down),]
hundred_genes_gost<-gost(query = hundred_genes$gene, organism = 'hsapiens',significant=T)
plot_hundred<-gostplot(hundred_genes_gost,capped = TRUE, interactive = T)
hundred_geens_enrich<-enrichGO(hundred_genes$gene_symbol, OrgDb='org.Hs.eg.db', keyType= "SYMBOL",ont= "MF",pAdjustMethod = "BH",pvalueCutoff  = 0.05,qvalueCutoff  = 0.05)
dotplot(hundred_geens_enrich)

###### what are the cytosolic genes ########

SP_cytosolic_down <- list_genesSP[list_genesSP$Symbol %in% cytosolic_genes_down_mine$gene,]
nrow(SP_cytosolic_down)

TM_cytosolic_down <- list_genesTM[list_genesTM$Symbol %in% cytosolic_genes_down_mine$gene,]
nrow(TM_cytosolic_down)

MITO_cytosolic_down <- result_df[result_df$Gene_Symbol %in% cytosolic_genes_down_mine$gene,]
nrow(MITO_cytosolic_down)

MITO_undefined_down <- result_df[result_df$Gene_Symbol %in% undefined_genes_down_mine$gene,]
nrow(MITO_undefined_down)

##### up genes ########

datax_8_up<-datax_8[datax_8$diffexpressed=="UP",]
genes_up_8<-enrichGO(datax_8_up$gene, OrgDb='org.Hs.eg.db', keyType= "SYMBOL",ont= "MF",pAdjustMethod = "BH",pvalueCutoff  = 0.05,qvalueCutoff  = 0.05)
dotplot(genes_up_8,showCategory=6)

up_cyto_8<-datax_8_up[datax_8_up$localisation=="cytosolic",]
genes_up_8<-enrichGO(up_cyto_8$gene, OrgDb='org.Hs.eg.db', keyType= "SYMBOL",ont= "MF",pAdjustMethod = "BH",pvalueCutoff  = 0.05,qvalueCutoff  = 0.05)
barplot(genes_up_8,showCategory =7)
edox <- setReadable(genes_up_8, 'org.Hs.eg.db', 'ENTREZID')
cnetplot(edox, node_label="all",color_category='firebrick', 
         color_gene='steelblue',cex_label_gene = 1,showCategory=7,categorySize="geneNum")


#datax_4 <- data.frame(LogFC_RNA = rna_results_4$log2FoldChange, Sig_RNA = as.numeric(rna_results_4$padj), gene = rna_results_4$gene)
#datax_4$Sig_RNA[is.na(datax_4$Sig_RNA)] <- 1
#temp_4 <- rep("NonSig", nrow(datax_4))
#temp_4[datax_4$Sig_RNA < 0.05 & abs(datax_4$LogFC_RNA) > 1] <- "Sig"
#datax_4$Sig <- factor(temp_4, levels = c("NonSig", "Sig"))

# add a column of NAs
#datax_4$diffexpressed <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
#datax_4$diffexpressed[datax_4$LogFC_RNA > 1 & datax_4$Sig == "Sig"] <- "UP"
#datax_4$diffexpressed[datax_4$LogFC_RNA < -1 & datax_4$Sig == "Sig"] <- "DOWN"

#datax_4_up<-datax_4[datax_4$diffexpressed=="UP",]
#genes_up_4<-enrichGO(datax_4_up$gene, OrgDb='org.Hs.eg.db', keyType= "SYMBOL",ont= "ALL",pAdjustMethod = "BH",pvalueCutoff  = 0.05,qvalueCutoff  = 0.05)
#dotplot(genes_up_4)

######### donw genes in general #########
#datax_4_down<-datax_4[datax_4$diffexpressed=="DOWN",]

membrane_genes_down_mine_4 <- datax_4_down[datax_4_down$gene %in% membrane_genes,]
cytosolic_genes_down_mine_4 <- datax_4_down[datax_4_down$gene %in% cytosolic_genes,]


genes_down_4<-enrichGO(datax_4_down$gene, OrgDb='org.Hs.eg.db', keyType= "SYMBOL",ont= "BP",pAdjustMethod = "BH",pvalueCutoff  = 0.05,qvalueCutoff  = 0.05)
dotplot(genes_down_4)

datax_8_downy<-datax_8_down[order(datax_8_down$LogFC_RNA,decreasing = F),]
ver<-datax_8_downy[1:500,]
datax_8_down<-datax_8[datax_8$diffexpressed=="DOWN",]
genes_down_8<-enrichGO(ver$gene, OrgDb='org.Hs.eg.db', keyType= "SYMBOL",ont= "CC",pAdjustMethod = "BH",pvalueCutoff  = 0.05,qvalueCutoff  = 0.05)
dotplot(genes_down_8)

#### 1. they dont contain SP or TM (TM only 1)

cytosolic_genes_ER<-gost(query = cytosolic_genes_down_mine$gene, organism = 'hsapiens',significant=T)
plot_cytosolic_ER<-gostplot(cytosolic_genes_ER,capped = TRUE, interactive = T)

cytosolic_ER_enrich<-enrichGO(cytosolic_genes_down_mine$gene, OrgDb='org.Hs.eg.db', keyType= "SYMBOL",ont= "MF",pAdjustMethod = "BH",pvalueCutoff  = 0.05,qvalueCutoff  = 0.05)
dotplot(cytosolic_ER_enrich,showCategory=5)
barplot(cytosolic_ER_enrich, showCategory=6) 
membrane_ER_enrich<-enrichGO(membrane_genes_down_mine$gene, OrgDb='org.Hs.eg.db', keyType= "SYMBOL",ont= "BP",pAdjustMethod = "BH",pvalueCutoff  = 0.05,qvalueCutoff  = 0.05)
dotplot(membrane_ER_enrich)

############ get protein lengths #############

combined_data_down_cyto_protein <- combined_data_plot[combined_data_plot$localisation=="cytosolic" & combined_data_plot$Set=="No SP/TM or mtRNA" & combined_data_plot$Source=="8 hours" & combined_data_plot$diffexpressed=="DOWN",]
nrow(combined_data_down_cyto_protein)
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
gene_names <- combined_data_down_cyto_protein$gene
# Use the listAttributes function to get valid attribute names
valid_attributes <- listAttributes(ensembl)

gene_info <- getBM(attributes = c("external_gene_name", "cds_length","transcript_length","transcript_biotype"), filters = "external_gene_name", values = gene_names, mart = ensembl)
protein_coding_genes <- gene_info[gene_info$transcript_biotype=="protein_coding",]

protein_long <- protein_coding_genes[protein_coding_genes$cds_length > 1000,]

# Extract unique gene names
unique_genes_RNA <- unique(protein_long$external_gene_name)
print(unique_genes_RNA)
# Count the number of unique genes
unique_gene_count_RNA <- length(unique_genes_RNA)

### From 193 proteins, 180 of those are considered long proteins, as they have more than 500 aminoacids

library(writexl)

# Define the file path for your Excel file
file_path <- "/Users/sgarcia/Desktop/protein_long.xlsx"

# Write the data frame to an Excel file
write_xlsx(protein_long, file_path)


######### what type are genes that variate ###########

lncRNA_down<-datax_8[datax_8$diffexpressed=="DOWN" & datax_8$gene_biotype=="lncRNA",]
lncRNA_up<-datax_8[datax_8$diffexpressed=="UP" & datax_8$gene_biotype=="lncRNA",]

############# what I am gonna do here is to take the down regulated samples at time 4 with no TM/SP
####### and see their structure

datax_4$diffexpressed[datax_4$LogFC_RNA < -1 & datax_4$Sig == "Sig"] <- "DOWN"

datax_4_down<-datax_4[datax_4$diffexpressed=="DOWN" & datax_4$Set=="No SP/TM or mtRNA",]

ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
gene_names_down <- datax_4_down$gene

gene_info <- getBM(attributes = c("external_gene_name", "cds_length","transcript_length","transcript_biotype"), filters = "external_gene_name", values = gene_names_down, mart = ensembl)
protein_coding_genes_down <- gene_info[gene_info$transcript_biotype=="protein_coding",]

protein_long_down <- protein_coding_genes_down[protein_coding_genes_down$cds_length > 1000,]

# Extract unique gene names
unique_genes_RNA <- unique(protein_long_down$external_gene_name)
print(unique_genes_RNA)
# Count the number of unique genes
unique_gene_count_RNA <- length(unique_genes_RNA)