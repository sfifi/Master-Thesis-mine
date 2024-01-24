
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

library(DESeq2)
library(EnsDb.Hsapiens.v86)
library(dplyr)
library(biomaRt)
library(ggplot2)



dir.path <- "/Volumes/AG_Landthaler/amanukyan/Projects/RiboSeq/Pipeline/run/sofia/output/FeatureCounts"
getwd()
# human data 
list_files <- list.files(dir.path)
list_files <- list_files[!grepl("summary", list_files)]

rnadata_cyto <- read.table(file="/Volumes/AG_Landthaler/amanukyan/Projects/RiboSeq/Pipeline/run/sofia_rna/output/count.txt", header = T, sep = "\t")
colnames(rnadata_cyto) <- lapply(colnames(rnadata_cyto), function(x) strsplit(x, split = "_S[0-9]+_")[[1]][1])
rownames(rnadata_cyto) <- rnadata_cyto[,1]
rnadata_cyto <- rnadata_cyto[,-1]
#rnadata_cyto <- rnadata_cyto[,c(6:13)]
rnadata_cyto <- rnadata_cyto[,c(6,7,9:13)]
colnames(rnadata_cyto) <- gsub("output.bam_files.", "", colnames(rnadata_cyto))




####
# DE analysis for cyto samples for three conditions (ctrl, 4h and 8 h) ####
####

# remove low count genes
count_gene_cyto <- apply(rnadata_cyto, 1, max) 
rnadata_cyto <- rnadata_cyto[which(count_gene_cyto > 10),]

rnaCond_cyto <- c("ctrl","ctrl","time 4 hours","time 4 hours","time 4 hours","time 8 hours","time 8 hours")
rnaCond_cyto<-factor(rnaCond_cyto)
# get riborex data
print(head(rnadata_cyto))


###### Compare ctrl 4 and 8 as done in ER ######

# get regular deseq2 results
dds <- DESeqDataSetFromMatrix(countData = rnadata_cyto,
                              DataFrame(rnaCond_cyto),
                              design= ~ rnaCond_cyto)
dds <- DESeq(dds)
unique_rnaCond_cyto <- as.character(rev(unique(rnaCond_cyto)))
rna_results_cyto <- NULL
for(i in 1:(length(unique_rnaCond_cyto)-1)){
  for(j in (i+1):(length(unique_rnaCond_cyto))){
    comparison <- c("rnaCond_cyto", unique_rnaCond_cyto[i], unique_rnaCond_cyto[j])
    print(comparison)
    cur_results <- as.data.frame(results(dds, contrast = comparison))
    rna_results_cyto <- rbind(rna_results_cyto, data.frame(cur_results, gene = rownames(cur_results), comparison = paste0(comparison[2:3], collapse = "_vs_")))
  }
}
rna_results_cyto$padj[is.na(rna_results_cyto$padj)] <- 1
rna_results_cyto_8<-filter(rna_results_cyto, comparison == "time 8 hours_vs_ctrl")
rna_results_cyto_4<-filter(rna_results_cyto, comparison == "time 4 hours_vs_ctrl")

# rna tpm
len_rna_cyto <- len[match(rownames(rnadata_cyto), len$hgnc_symbol),]
len_rna_cyto <- na.omit(len_rna_cyto)
rnadata_tpm_cyto <- rnadata_cyto[intersect(rownames(rnadata_cyto), len_rna_cyto$hgnc_symbol),]
rnadata_tpm_cyto <- sweep(rnadata_tpm_cyto * 1000, 1, len_rna_cyto$len, FUN = "/")
rnadata_tpm_cyto <- sweep(rnadata_tpm_cyto * 1000000, 2, colSums(rnadata_tpm_cyto), FUN = "/")

# calculate log RNA mean for each condition
cond_names <- rnaCond_cyto
Mean_RNA_Data_cyto <- aggregate(t(rnadata_tpm_cyto), list(cond_names), mean)
Mean_RNA_Data_cyto <- t(Mean_RNA_Data_cyto)
colnames(Mean_RNA_Data_cyto) <- Mean_RNA_Data_cyto[1,]
Mean_RNA_Data_cyto <- data.frame(apply(Mean_RNA_Data_cyto[-1,],2,as.numeric), row.names = rownames(Mean_RNA_Data_cyto[-1,]))
LogMean_RNA_Data_cyto <- log(Mean_RNA_Data_cyto + 0.001)

datax_cyto_4 <- data.frame(LogFC_RNA = rna_results_cyto_4$log2FoldChange, Sig_RNA = as.numeric(rna_results_cyto_4$padj), gene = rna_results_cyto_4$gene)
datax_cyto_4$Sig_RNA[is.na(datax_cyto_4$Sig_RNA)] <- 1
temp_cyto_4 <- rep("NonSig", nrow(datax_cyto_4))
temp_cyto_4[datax_cyto_4$Sig_RNA < 0.05 & abs(datax_cyto_4$LogFC_RNA) > 0.5] <- "Sig"
datax_cyto_4$Sig <- factor(temp_cyto_4, levels = c("NonSig", "Sig"))
cur_LogMean_RNA_Data_cyto_4 <- LogMean_RNA_Data_cyto[datax_cyto_4$gene,]

cur_LogMean_RNA_Data_cyto_4_select <- data.frame(x = as.numeric(cur_LogMean_RNA_Data_cyto_4[,1]), 
                                                  y = as.numeric(cur_LogMean_RNA_Data_cyto_4[,2]), 
                                                  Sig = as.character(datax_cyto_4$Sig))

datax_cyto_8 <- data.frame(LogFC_RNA = rna_results_cyto_8$log2FoldChange, Sig_RNA = as.numeric(rna_results_cyto_8$padj), gene = rna_results_cyto_8$gene)
datax_cyto_8$Sig_RNA[is.na(datax_cyto_8$Sig_RNA)] <- 1
temp_cyto_8 <- rep("NonSig", nrow(datax_cyto_8))
temp_cyto_8[datax_cyto_8$Sig_RNA < 0.05 & abs(datax_cyto_8$LogFC_RNA) > 0.5] <- "Sig"
datax_cyto_8$Sig <- factor(temp_cyto_8, levels = c("NonSig", "Sig"))
cur_LogMean_RNA_Data_cyto_8 <- LogMean_RNA_Data_cyto[datax_cyto_8$gene,]

cur_LogMean_RNA_Data_cyto_8_select <- data.frame(x = as.numeric(cur_LogMean_RNA_Data_cyto_8[,1]), 
                                                  y = as.numeric(cur_LogMean_RNA_Data_cyto_8[,3]), 
                                                  Sig = as.character(datax_cyto_8$Sig))

# add a column of NAs
cur_LogMean_RNA_Data_cyto_4_select$diffexpressed <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
cur_LogMean_RNA_Data_cyto_4_select$diffexpressed[datax_cyto_4$LogFC_RNA > 0.5 & datax_cyto_4$Sig == "Sig"] <- "UP"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
cur_LogMean_RNA_Data_cyto_4_select$diffexpressed[datax_cyto_4$LogFC_RNA < -0.5 & datax_cyto_4$Sig == "Sig"] <- "DOWN"

cur_LogMean_RNA_Data_cyto_4_select$gene_symbol<-datax_cyto_4$gene

genes_up_RNA_cyto_4<-cur_LogMean_RNA_Data_cyto_4_select%>%
  dplyr::filter(diffexpressed == "UP")%>%
  dplyr::filter(Sig=="Sig")%>%
  select(gene_symbol)

membrane_cyto<-genes_up_RNA_cyto_4[genes_up_RNA_cyto_4$gene_symbol %in% membrane_genes,]
membrane_cyto<-sample(membrane_cyto,10)
membrane_cyto<-c(membrane_cyto,"AGO2")
#labelling
cur_LogMean_RNA_Data_cyto_4_select$delabel <- NA
cur_LogMean_RNA_Data_cyto_4_select$delabel[cur_LogMean_RNA_Data_cyto_4_select$gene_symbol %in% membrane_cyto] <- cur_LogMean_RNA_Data_cyto_4_select$gene_symbol[cur_LogMean_RNA_Data_cyto_4_select$gene_symbol %in% membrane_cyto]


# add a column of NAs
cur_LogMean_RNA_Data_cyto_8_select$diffexpressed <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
cur_LogMean_RNA_Data_cyto_8_select$diffexpressed[datax_cyto_8$LogFC_RNA > 0.5 & datax_cyto_8$Sig == "Sig"] <- "UP"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
cur_LogMean_RNA_Data_cyto_8_select$diffexpressed[datax_cyto_8$LogFC_RNA < -0.5 & datax_cyto_8$Sig == "Sig"] <- "DOWN"

cur_LogMean_RNA_Data_cyto_8_select$gene_symbol<-datax_cyto_8$gene

genes_up_RNA_cyto_8<-cur_LogMean_RNA_Data_cyto_8_select%>%
  dplyr::filter(diffexpressed == "UP")%>%
  dplyr::filter(Sig=="Sig")%>%
  select(gene_symbol)
membrane_cyto8<-genes_up_RNA_cyto_8[genes_up_RNA_cyto_8$gene_symbol %in% membrane_genes,]
membrane_cyto8<-sample(membrane_cyto8,10)
membrane_cyto8<-c(membrane_cyto8,"AGO2")

#labelling
cur_LogMean_RNA_Data_cyto_8_select$delabel <- NA
cur_LogMean_RNA_Data_cyto_8_select$delabel[cur_LogMean_RNA_Data_cyto_8_select$gene_symbol %in% membrane_cyto8] <- cur_LogMean_RNA_Data_cyto_8_select$gene_symbol[cur_LogMean_RNA_Data_cyto_8_select$gene_symbol %in% membrane_cyto8]



genes_down_RNA_cyto_8<-cur_LogMean_RNA_Data_cyto_8_select%>%
  dplyr::filter(diffexpressed == "DOWN")%>%
  select(gene_symbol)

mycolors <- c("blue", "red", "grey")
names(mycolors) <- c("DOWN", "UP", "NO")

#### Visulaizationn of downregulated genes after 4 and 8 hours of auxin depletion


p5<-ggplot(cur_LogMean_RNA_Data_cyto_4_select, aes(x = x, y = y,label=delabel,col=diffexpressed)) + 
  geom_point() +   xlab("Non treated") +
  ylab("4 hours Auxin treatment") + geom_text_repel(max.overlaps=40,size=5,nudge_y = 3) + scale_colour_manual(values = mycolors) + theme_minimal() +
  theme(
    axis.title.x = element_text(size = rel(2)), # make x axis label 1.5 times larger
    axis.title.y = element_text(size = rel(2)), # make y axis label 1.5 times larger
    plot.title = element_text(size = rel(2)),
    legend.text = element_text(size = rel(1)),
    legend.title = element_text(size = rel(1.5))# make title 2 times larger
  )

p6<-ggplot(cur_LogMean_RNA_Data_cyto_8_select, aes(x = x, y = y,label=delabel,col=diffexpressed)) + 
  geom_point() +   xlab("Non treated") +
  ylab("8 hours Auxin treatment") + geom_text_repel(max.overlaps = 40,size=5,nudge_y = 3) + scale_colour_manual(values = mycolors) + theme_minimal() +
  theme(
    axis.title.x = element_text(size = rel(2)), # make x axis label 1.5 times larger
    axis.title.y = element_text(size = rel(2)), # make y axis label 1.5 times larger
    plot.title = element_text(size = rel(2)),
    legend.text = element_text(size = rel(1)),
    legend.title = element_text(size = rel(1.5))# make title 2 times larger
  )
#arrangment in a grid of the previous plots

layout <- rbind(c(1, 2))
grid.arrange(p5, p6,layout_matrix=layout)

############## Sig

# add a column of NAs
datax_cyto_4$diffexpressed <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
datax_cyto_4$diffexpressed[datax_cyto_4$LogFC_RNA > 1 & datax_cyto_4$Sig == "Sig"] <- "UP"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
datax_cyto_4$diffexpressed[datax_cyto_4$LogFC_RNA < -1 & datax_cyto_4$Sig == "Sig"] <- "DOWN"


# add a column of NAs
datax_cyto_8$diffexpressed <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
datax_cyto_8$diffexpressed[datax_cyto_8$LogFC_RNA > 1 & datax_cyto_8$Sig == "Sig"] <- "UP"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
datax_cyto_8$diffexpressed[datax_cyto_8$LogFC_RNA < -1 & datax_cyto_8$Sig == "Sig"] <- "DOWN"


############# SET separation ####################
genes_MITO<-datax_cyto_4$gene %in% result_df$Gene_Symbol
genes_SP<-datax_cyto_4$gene %in% list_genesSP$Symbol
genes_TM<-datax_cyto_4$gene %in% list_genesTM$Symbol
datax_cyto_4$Set <- NA
datax_cyto_4$Set[genes_MITO] <- "mtRNA"
datax_cyto_4$Set[genes_SP] <- "SP"
datax_cyto_4$Set[genes_TM] <- "TM"
# Define the categories you want to exclude
categories_to_exclude <- c("mtRNA", "SP", "TM")
# Use ifelse to update the "Set" column
datax_cyto_4$Set <- ifelse(datax_cyto_4$Set %in% categories_to_exclude, datax_cyto_4$Set, "No SP/TM or mtRNA")


genes_MITO<-datax_cyto_8$gene %in% result_df$Gene_Symbol
genes_SP<-datax_cyto_8$gene %in% list_genesSP$Symbol
genes_TM<-datax_cyto_8$gene %in% list_genesTM$Symbol
datax_cyto_8$Set <- NA
datax_cyto_8$Set[genes_MITO] <- "mtRNA"
datax_cyto_8$Set[genes_SP] <- "SP"
datax_cyto_8$Set[genes_TM] <- "TM"
# Define the categories you want to exclude
categories_to_exclude <- c("mtRNA", "SP", "TM")
# Use ifelse to update the "Set" column
datax_cyto_8$Set <- ifelse(datax_cyto_8$Set %in% categories_to_exclude, datax_cyto_8$Set, "No SP/TM or mtRNA")

###### localisation #########

datax_cyto_8$localisation <- NA
datax_cyto_8$localisation[datax_cyto_8$gene %in% membrane_genes] <- "membrane"
datax_cyto_8$localisation[datax_cyto_8$gene %in% cytosolic_genes] <- "cytosolic"
condition <- datax_cyto_8$gene %in% undefined_genes | is.na(datax_cyto_8$localisation)
datax_cyto_8$localisation[condition] <- "undefined"

datax_cyto_4$localisation <- NA
datax_cyto_4$localisation[datax_cyto_4$gene %in% membrane_genes] <- "membrane"
datax_cyto_4$localisation[datax_cyto_4$gene %in% cytosolic_genes] <- "cytosolic"
condition <- datax_cyto_4$gene %in% undefined_genes | is.na(datax_cyto_4$localisation)
datax_cyto_4$localisation[condition] <- "undefined"


genes_MITO<-datax_cyto_4$gene %in% result_df$Gene_Symbol
genes_SP<-datax_cyto_4$gene %in% list_genesSP$Symbol
genes_TM<-datax_cyto_4$gene %in% list_genesTM$Symbol
datax_cyto_4$Set <- NA
datax_cyto_4$Set[genes_MITO] <- "mtRNA"
datax_cyto_4$Set[genes_SP] <- "SP"
datax_cyto_4$Set[genes_TM] <- "TM"
categories_to_exclude <- c("mtRNA", "SP", "TM")
datax_cyto_4$Set <- ifelse(datax_cyto_4$Set %in% categories_to_exclude, datax_cyto_4$Set, "No SP/TM or mtRNA")

genes_MITO<-datax_cyto_8$gene %in% result_df$Gene_Symbol
genes_SP<-datax_cyto_8$gene %in% list_genesSP$Symbol
genes_TM<-datax_cyto_8$gene %in% list_genesTM$Symbol
datax_cyto_8$Set <- NA
datax_cyto_8$Set[genes_MITO] <- "mtRNA"
datax_cyto_8$Set[genes_SP] <- "SP"
datax_cyto_8$Set[genes_TM] <- "TM"
# Define the categories you want to exclude
categories_to_exclude <- c("mtRNA", "SP", "TM")
# Use ifelse to update the "Set" column
datax_cyto_8$Set <- ifelse(datax_cyto_8$Set %in% categories_to_exclude, datax_cyto_8$Set, "No SP/TM or mtRNA")

###### addition of byotype ###########

genes_cyto_4 <- datax_cyto_4$gene 
gene_info_cyto_4 <- getBM(
  attributes = c("external_gene_name", "gene_biotype"),
  filters = "external_gene_name",
  values = genes_cyto_4,
  mart = ensembl
)
datax_cyto_4 <- datax_cyto_4 %>%
  left_join(gene_info_cyto_4, by = c("gene" = "external_gene_name"))%>%
  mutate(gene_biotype = coalesce(gene_biotype, "noinfo"))

genes_cyto_8 <- datax_cyto_8$gene 
gene_info_cyto_8 <- getBM(
  attributes = c("external_gene_name", "gene_biotype"),
  filters = "external_gene_name",
  values = genes_cyto_8,
  mart = ensembl
)
datax_cyto_8 <- datax_cyto_8 %>%
  left_join(gene_info_cyto_8, by = c("gene" = "external_gene_name"))%>%
  mutate(gene_biotype = coalesce(gene_biotype, "noinfo"))

######## Source addition

# Add a column to datax_8 to indicate the source
datax_cyto_8 <- datax_cyto_8 %>%
  mutate(Source = "datax_8")
# Add a column to datax_4 to indicate the source
datax_cyto_4 <- datax_cyto_4 %>%
  mutate(Source = "datax_4")


# Step 1: Filter significant genes from datax_8
significant_datax_8 <- subset(datax_cyto_8, Sig == "Sig")
# Step 2: Get corresponding genes from datax_4
corresponding_datax_4 <- datax_cyto_4[datax_cyto_4$gene %in% significant_datax_8$gene,]
# Step 3: Combine the two sets of genes
genes_for_plotting <- rbind(significant_datax_8, corresponding_datax_4)
# Step 1: Filter significant genes from datax_4
significant_datax_4 <- subset(datax_cyto_4, Sig == "Sig")
# Step 2: Get corresponding genes from datax_8
corresponding_datax_8 <- datax_cyto_8[datax_cyto_8$gene %in% significant_datax_4$gene,]
# Step 3: Combine the two sets of genes
genes_for_plotting_reversed <- rbind(significant_datax_4, corresponding_datax_8)
combined_genes_for_plotting <- rbind(genes_for_plotting, genes_for_plotting_reversed)
combined_genes_for_plotting<- combined_genes_for_plotting %>%
  mutate(Source = ifelse(Source == "datax_4", "4 h", ifelse(Source == "datax_8", "8 h", Source)))

combined_data_plot <-  combined_genes_for_plotting
combined_data_plot$Group <- paste0(combined_data_plot$Set, ' (', combined_data_plot$Source, ')')
combined_data_plot <- combined_data_plot %>%
  filter(localisation %in% c("cytosolic", "membrane") & localisation != "undefined")
count_data <- combined_data_plot %>%
  group_by(Group,localisation) %>%
  summarise(Count = n(), .groups = 'drop')


# Create the violin plot UP

ggplot(combined_data_plot, aes(x = Group, y = LogFC_RNA)) +
  geom_violin(aes(fill = Set),alpha=0.5) + geom_line(aes(group=gene), linewidth=0.4, color="gray") +
  geom_hline(yintercept = 0) + 
  geom_jitter(data=combined_data_plot, aes(x=Group, y=LogFC_RNA), width=.1, alpha=.3) +
  labs(title = "LogFC RNA on cytosol", x = "Sets", y = "LogFC RNA") +
  scale_fill_manual(values = c("SP" = "orange", "TM" = "dark green", "mtRNA" = "purple", "No SP/TM or mtRNA" = "grey")) +
  theme_minimal() +
  geom_text(data = count_data, aes(x = Group, y = Inf, label = Count), vjust = 1.5, size = 3, position = position_jitter(width = 0.1)) + 
  ylim(-5,5) + facet_grid(~localisation) + theme(
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


########## simplified plot ###########

# Filter the data to include only "cytosolic" and "membrane" localizations
combined_data_plot <-  combined_genes_for_plotting
filtered_data_cyto <- combined_data_plot %>%
  filter(localisation %in% c("cytosolic", "membrane") & localisation != "undefined")
count_data <-filtered_data_cyto %>%
  group_by(localisation,Source) %>%
  summarise(Count = n(), .groups = 'drop')


# Create the plot with the filtered data without using "Group"
ggplot(filtered_data_cyto, aes(x = Source, y = LogFC_RNA)) +
  geom_violin(aes(fill = localisation), alpha = 0.5) + 
  geom_line(aes(group = gene), linewidth = 0.4, color = "gray") +
  geom_hline(yintercept = 0) + 
  geom_jitter(data = filtered_data_cyto, aes(x = Source, y = LogFC_RNA), width = .1, alpha = .3) +
  labs(title = "LogFC RNA on cytosol", x = "Time points", y = "LogFC RNA") +
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



######### whats up whats down #############
#########################################
#######################################

###### UP
genes_up_RNA_cyto_8$localisation <- NA
genes_up_RNA_cyto_8$localisation[genes_up_RNA_cyto_8$gene %in% membrane_genes] <- "membrane"
genes_up_RNA_cyto_8$localisation[genes_up_RNA_cyto_8$gene %in% cytosolic_genes] <- "cytosolic"
genes_up_RNA_cyto_8$localisation[genes_up_RNA_cyto_8$gene %in% undefined_genes] <- "undefined"
genes_up_RNA_cyto_8<-genes_up_RNA_cyto_8[genes_up_RNA_cyto_8$localisation=="membrane",]

genes_up_RNA_cyto_4$localisation <- NA
genes_up_RNA_cyto_4$localisation[genes_up_RNA_cyto_4$gene %in% membrane_genes] <- "membrane"
genes_up_RNA_cyto_4$localisation[genes_up_RNA_cyto_4$gene %in% cytosolic_genes] <- "cytosolic"
genes_up_RNA_cyto_4$localisation[genes_up_RNA_cyto_4$gene %in% undefined_genes] <- "undefined"
genes_up_RNA_cyto_4<-genes_up_RNA_cyto_4[genes_up_RNA_cyto_4$localisation=="cytosolic",]


#dotplot
up_cyto_8_enrichGO<-enrichGO(genes_up_RNA_cyto_8$gene_symbol, OrgDb='org.Hs.eg.db', keyType= "SYMBOL",ont= "BP",pAdjustMethod = "BH",pvalueCutoff  = 0.05,qvalueCutoff  = 0.05)
dotplot(up_cyto_8_enrichGO)
#gost
up_gost_cyto_8<-gost(query = genes_up_RNA_cyto_8, organism = 'hsapiens',significant=T)
up_cyto_gostplot<-gostplot(up_gost_cyto_8,capped = TRUE, interactive = T)
#gost
up_cyto<-datax_cyto_4[datax_cyto_4$diffexpressed=="UP" & datax_cyto_4$Sig=="Sig",]
data <- up_cyto %>%
  arrange(desc(LogFC_RNA))
data_select<-data[1:800,]
up_cyto_4_enrichGO<-enrichGO(data_select$gene, OrgDb='org.Hs.eg.db', keyType= "SYMBOL",ont= "MF",pAdjustMethod = "BH",pvalueCutoff  = 0.05,qvalueCutoff  = 0.05)
barplot(up_cyto_4_enrichGO,showCategory=8)
cnetplot(up_cyto_4_enrichGO, node_label="all",color_category='firebrick', 
         color_gene='steelblue',cex_label_gene = 1,showCategory=8,categorySize="geneNum")
####### DOWN

genes_down_RNA_cyto_4$localisation <- NA
genes_down_RNA_cyto_4$localisation[genes_down_RNA_cyto_4$gene %in% membrane_genes] <- "membrane"
genes_down_RNA_cyto_4$localisation[genes_down_RNA_cyto_4$gene %in% cytosolic_genes] <- "cytosolic"
genes_down_RNA_cyto_4$localisation[genes_down_RNA_cyto_4$gene %in% undefined_genes] <- "undefined"
genes_down_RNA_cyto_4<-genes_down_RNA_cyto_4[genes_down_RNA_cyto_4$localisation=="cytosolic",]
#dotplot
down_cyto<-datax_cyto_4[datax_cyto_4$diffexpressed=="DOWN" & datax_cyto_4$Sig=="Sig",]
data <- down_cyto %>%
  arrange(LogFC_RNA)
data_select<-data[1:300,]

down_cyto_8_enrichGO<-enrichGO(data_select$gene, OrgDb='org.Hs.eg.db', keyType= "SYMBOL",ont= "CC",pAdjustMethod = "BH",pvalueCutoff  = 0.05,qvalueCutoff  = 0.05)
barplot(down_cyto_8_enrichGO,showCategory=4)
down_cyto_4_enrichGO<-enrichGO(genes_down_RNA_cyto_4$gene_symbol, OrgDb='org.Hs.eg.db', keyType= "SYMBOL",ont= "BP",pAdjustMethod = "BH",pvalueCutoff  = 0.05,qvalueCutoff  = 0.05)
dot2<-barplot(down_cyto_4_enrichGO,showCategory=4)
edox <- setReadable(down_cyto_8_enrichGO, 'org.Hs.eg.db', 'ENTREZID')
cnetplot(down_cyto_8_enrichGO, node_label="all",color_category='firebrick', 
         color_gene='steelblue',cex_label_gene = 1,showCategory=4,categorySize="geneNum")

#puttogether
layout <- rbind(c(1, 2))
grid.arrange(dot2, dot1,layout_matrix=layout)
#gost
down_gost_cyto_8<-gost(query = genes_down_RNA_cyto_8$gene_symbol, organism = 'hsapiens',significant=T)
down_cyto_gostplot<-gostplot(down_gost_cyto_8,capped = TRUE, interactive = T)


######### SP ############
##########################
#######################

list_genesER<-read.table(file="/Volumes/AG_Landthaler/sofia/riboSRP72/sp_list.txt", header = T, sep = "\t")
genes_SP_up_cyto<-list_genesER[list_genesER$Symbol %in% genes_up_RNA_cyto_8$gene_symbol,]
genes_SP_up_cyto<-print(genes_SP_up_cyto$Symbol)
sp_cyto_up<-enrichGO(genes_SP_up_cyto, OrgDb='org.Hs.eg.db', keyType= "SYMBOL",ont= "ALL",pAdjustMethod = "BH",pvalueCutoff  = 0.05,qvalueCutoff  = 0.05)
dotplot(sp_cyto_up)

cyto_up_SP<-gost(query = genes_SP_up_cyto, organism = 'hsapiens',significant=T)
cyto_enrich_SP<-gostplot(cyto_up_SP,capped = TRUE, interactive = T)

############ TM ################
##############################
#############################

list_genesTM<-read.table(file="/Users/sgarcia/Downloads/mmilek-hdlbp_rev-c0e27b8/data/tm_list.txt", header = T, sep = "\t")
genes_TM_up_8_cyto<-list_genesTM[list_genesTM$Symbol %in% genes_up_RNA_cyto_8$gene_symbol,]
genes_TM_up_8_cyto<-print(genes_TM_up_8_cyto$Symbol)

tm_cyto_up<-enrichGO(genes_TM_up_cyto, OrgDb='org.Hs.eg.db', keyType= "SYMBOL",ont= "ALL",pAdjustMethod = "BH",pvalueCutoff  = 0.05,qvalueCutoff  = 0.05)
dotplot(tm_cyto_up)

cyto_up_TM<-gost(query = genes_TM_up_cyto, organism = 'hsapiens',significant=T)
cyto_enrich_TM<-gostplot(cyto_up_TM,capped = TRUE, interactive = T)

######### MITO: run before genes_cyto_RNA_analyses.R ##########

genes_MITO_down_8_cyto<-result_df[result_df$Gene_Symbol %in% genes_down_RNA_cyto_4$gene_symbol,]



######## comparison ER vs cyto ############

genes_down_RNA_cyto_4<-cur_LogMean_RNA_Data_cyto_4_select%>%
  dplyr::filter(diffexpressed == "DOWN")%>%
  select(gene_symbol)


rnadata_comp<-rnadata_tpm_cyto[rnadata_tpm_cyto$gene %in% datax_8_down$gene,]

comparison_ER_total_8_down <- genes_down_RNA_cyto_8[genes_down_RNA_cyto_8$gene_symbol %in% datax_8_down$gene,]
comparison_ER_total_8_up <- as.data.frame(genes_up_RNA_cyto_8[genes_up_RNA_cyto_8$gene_symbol %in% datax_8_down$gene,])

comparison_ER_total_4_down<- genes_down_RNA_cyto_4[genes_down_RNA_cyto_4$gene_symbol %in% datax_4_down$gene,]
comparison_ER_total_4_up<- genes_up_RNA_cyto_4[genes_up_RNA_cyto_4$gene_symbol %in% datax_4_down$gene,]

membrane_in_cytosol<-genes_up_RNA_cyto_8[genes_up_RNA_cyto_8$gene %in% membrane_genes$gene,]
membrane_in_cytosol_4<-genes_up_RNA_cyto_4[genes_up_RNA_cyto_4$gene_symbol %in% membrane_genes,]
cytosolic_in_cytosol_4<-genes_up_RNA_cyto_4[genes_up_RNA_cyto_4$gene_symbol %in% cytosolic_genes,]

membrane<-datax_cyto_8[datax_cyto_8$gene %in% membrane_genes$gene,]

rnadata_tpm_cyto$gene<-rownames(rnadata_tpm_cyto)


rnadata_comp_cyto<-datax_cyto_8[datax_cyto_8$gene %in% datax_8_down$gene,]
see<-rnadata_comp_cyto[rnadata_comp_cyto$Sig_RNA<0.9 & rnadata_comp_cyto$LogFC_RNA>0.5,]
obs<-datax_cyto_4[datax_cyto_4$gene %in% comparison_ER_total_4_up,]


no<-datax_cyto_8[datax_cyto_8$Set=="No SP/TM or mtRNA" & datax_cyto_8$Sig=="Sig" & datax_cyto_8$diffexpressed=="UP",]
