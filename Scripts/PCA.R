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

rnadata_ER <- read.table(file="/Volumes/AG_Landthaler/amanukyan/Projects/RiboSeq/Pipeline/run/sofia_rna/output/count.txt", header = T, sep = "\t")
colnames(rnadata_ER) <- lapply(colnames(rnadata_ER), function(x) strsplit(x, split = "_S[0-9]+_")[[1]][1])
rownames(rnadata_ER) <- rnadata_ER[,1]
rnadata_ER <- rnadata_ER[,-1]
rnadata_ER <- rnadata_ER[,c(14:21)]
colnames(rnadata_ER) <- gsub("output.bam_files.", "", colnames(rnadata_ER))

rnadata_cyto <- read.table(file="/Volumes/AG_Landthaler/amanukyan/Projects/RiboSeq/Pipeline/run/sofia_rna/output/count.txt", header = T, sep = "\t")
colnames(rnadata_cyto) <- lapply(colnames(rnadata_cyto), function(x) strsplit(x, split = "_S[0-9]+_")[[1]][1])
rownames(rnadata_cyto) <- rnadata_cyto[,1]
rnadata_cyto <- rnadata_cyto[,-1]
rnadata_cyto <- rnadata_cyto[,c(6:13)]
colnames(rnadata_cyto) <- gsub("output.bam_files.", "", colnames(rnadata_cyto))

total_fractions<- read.table(file="/Volumes/AG_Landthaler/amanukyan/Projects/RiboSeq/Pipeline/run/sofia_rna/output/count.txt", header = T, sep = "\t")
colnames(total_fractions) <- lapply(colnames(total_fractions), function(x) strsplit(x, split = "_S[0-9]+_")[[1]][1])
rownames(total_fractions) <- total_fractions[,1]
total_fractions <- total_fractions[,c(23:30)]
colnames(total_fractions) <- gsub("output.bam_files.", "", colnames(total_fractions))
colnames(total_fractions) <- gsub(".bam", "", colnames(total_fractions))

#remove low count genes
count_gene <- apply(rnadata_ER, 1, max) 
rnadata_ER <- rnadata_ER[which(count_gene > 10),]

count_gene <- apply(rnadata_cyto, 1, max) 
rnadata_cyto <- rnadata_cyto[which(count_gene > 10),]

count_gene <- apply(total_fractions, 1, max) 
total_fractions <- total_fractions[which(count_gene > 10),]


rnaCond <- c("ctrl","ctrl","ctrl","time 4 hours","time 4 hours","time 4 hours","time 8 hours","time 8 hours")

rnaCond<-factor(rnaCond)


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
len_rna <- len[match(rownames(rnadata_ER), len$hgnc_symbol),]
len_rna <- na.omit(len_rna)
rnadata_tpm_ER <- rnadata_ER[intersect(rownames(rnadata_ER), len_rna$hgnc_symbol),]
rnadata_tpm_ER <- sweep(rnadata_tpm_ER * 1000, 1, len_rna$len, FUN = "/")
rnadata_tpm_ER <- sweep(rnadata_tpm_ER * 1000000, 2, colSums(rnadata_tpm_ER), FUN = "/")

len_rna <- len[match(rownames(rnadata_cyto), len$hgnc_symbol),]
len_rna <- na.omit(len_rna)
rnadata_tpm_cyto <- rnadata_cyto[intersect(rownames(rnadata_cyto), len_rna$hgnc_symbol),]
rnadata_tpm_cyto <- sweep(rnadata_tpm_cyto * 1000, 1, len_rna$len, FUN = "/")
rnadata_tpm_cyto <- sweep(rnadata_tpm_cyto * 1000000, 2, colSums(rnadata_tpm_cyto), FUN = "/")

len_rna <- len[match(rownames(total_fractions), len$hgnc_symbol),]
len_rna <- na.omit(len_rna)
rnadata_tpm <- total_fractions[intersect(rownames(total_fractions), len_rna$hgnc_symbol),]
rnadata_tpm <- sweep(rnadata_tpm * 1000, 1, len_rna$len, FUN = "/")
rnadata_tpm <- sweep(rnadata_tpm * 1000000, 2, colSums(rnadata_tpm), FUN = "/")


library(debrowser)

normdata_ER <- getNormalizedMatrix(rnadata_tpm_ER, method = "TMM")

normdata_cyto <- getNormalizedMatrix(rnadata_tpm_cyto, method = "TMM")

normdata<- getNormalizedMatrix(rnadata_tpm, method = "TMM")

normdata_pr_ER <- apply(normdata_ER,1,scale)
prtemp_ER <- prcomp(normdata_pr_ER)
normdata_pr_ER <- prtemp_ER$x
# inclusion to see the names of the columns
col<-colnames(normdata_pr_ER)
ggplot(data.frame(normdata_pr_ER), aes(x = PC1, y = PC2, color = rnaCond)) + geom_point() + labs(title = "PCA") + geom_text(aes(label = col, vjust = 1.5, hjust = 0.5))


normdata_pr_cyto <- apply(normdata_cyto,1,scale)
prtemp_cyto <- prcomp(normdata_pr_cyto)
normdata_pr_cyto <- prtemp_cyto$x
colc<-colnames(normdata_pr_cyto)
ggplot(data.frame(normdata_pr_cyto), aes(x = PC1, y = PC2, color = rnaCond)) + geom_point() + labs(title = "PCA") + geom_text(aes(label = colc, vjust = 1.5, hjust = 0.5))

normdata_pr <- apply(normdata,1,scale)
prtemp <- prcomp(normdata_pr)
normdata_pr <- prtemp$x
colt<-colnames(normdata_pr)
ggplot(data.frame(normdata_pr), aes(x = PC1, y = PC2, color = rnaCond)) + geom_point() + labs(title = "PCA") + geom_text(aes(label = colt, vjust = 1.5, hjust = 0.5))

####
# DE analysis ####
####

#Differential expression of RNA on ER among conditions
dds <- DESeqDataSetFromMatrix(countData = rnadata_ER,
                              DataFrame(rnaCond),
                              design= ~ rnaCond)
dds <- DESeq(dds)
unique_rnaCond <- as.character(rev(unique(rnaCond)))
rna_results <- NULL
for(i in 1:(length(unique_riboCond)-1)){
  for(j in (i+1):(length(unique_rnaCond))){
    comparison <- c("rnaCond", unique_rnaCond[i], unique_rnaCond[j])
    print(comparison)
    cur_results <- as.data.frame(results(dds, contrast = comparison))
    rna_results <- rbind(rna_results, data.frame(cur_results, gene = rownames(cur_results), comparison = paste0(comparison[2:3], collapse = "_vs_")))
  }
}
rna_results$padj[is.na(rna_results$padj)] <- 1


#For each gene, we count the total number of reads for that gene in all samples 
#and remove those that don't have at least 1 read. 
dds <- dds[ rowSums(DESeq2::counts(dds)) > 1, ]
dds <- DESeq(dds)

#compute the contrast for the 'group' variable where 'CTRL' 
#samples are used as the control group. 
DEresults <- results(dds, contrast = comparison)
#sort results by increasing p-value

DEresults$padj[is.na(DEresults$padj)] <- 1

DEresults <- DEresults[DEresults$padj<0.5,]
DEresults_1 <- DEresults[DEresults$log2FoldChange<1.5,]
genes<-rownames(DEresults_1)



#shows a summary of the results
print(DEresults)

DESeq2::plotMA(object = dds, ylim = c(-5, 5))

ggplot(data = as.data.frame(DEresults), aes(x = pvalue)) + 
  geom_histogram(bins = 100)

# extract normalized counts from the DESeqDataSet object
countsNormalized <- DESeq2::counts(dds, normalized = TRUE)
countsnoNormalized <- DESeq2::counts(dds, normalized = FALSE)

rld <- rlog(dds)
DESeq2::plotPCA(rld, ntop = 500, intgroup = "rnaCond") 

################### TOTAL COUNTS ASSIGNED BAM FILES ######################


reads <- read.table(file="/Volumes/AG_Landthaler/sofia/riboSRP72/Ribo/FeatureCounts/cds_counts.txt.summary", header = T, sep = "\t")
colnames(reads) <- gsub("output.bam_files.", "", colnames(reads))
colnames(reads) <- gsub(".bam$", "", colnames(reads))

reads_assigned <- reads[1,]

reads_assigned <- reads_assigned[,-1]

reads_assigned <- t(reads_assigned)

colnames(reads_assigned) <- "Status"

reads_assigned <- as.data.frame(reads_assigned)

reads_assigned$Samples <- rownames(reads_assigned)


# Create a barplot using ggplot2
ggplot(reads_assigned, aes(x = Samples, y = Status, fill = substr(Samples, 1, 1))) +
  geom_bar(stat = "identity") +
  labs(
    title = "CDS ribosome profiling counts assigned to BAM files",
    x = "Samples",
    y = "Counts"
  ) +
  scale_fill_manual(
    values = c("T" = "blue", "C" = "orange", "E" = "brown"),
    name = "Localization",labels = c("Total" = "Tissue", "C" = "cytosol", "E" = "ER")) + 
  
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

reads_rna <- read.table(file="/Volumes/AG_Landthaler/sofia/riboSRP72/RNA/count.txt.summary", header = T, sep = "\t")
colnames(reads_rna) <- gsub("output.bam_files.SRP72_", "", colnames(reads_rna))
colnames(reads_rna) <- gsub("_sorted.bam$", "", colnames(reads_rna))

reads_assigned_rna <- reads_rna[1,]

reads_assigned_rna <- reads_assigned_rna[,-1]

reads_assigned_rna <- t(reads_assigned_rna)

colnames(reads_assigned_rna) <- "Status"

reads_assigned_rna <- as.data.frame(reads_assigned_rna)

reads_assigned_rna$Samples <- rownames(reads_assigned_rna)


# Create a barplot using ggplot2
ggplot(reads_assigned_rna, aes(x = Samples, y = Status, fill = substr(Samples, 1, 1))) +
  geom_bar(stat = "identity") +
  labs(
    title = "mRNA counts assigned to BAM files",
    x = "Samples",
    y = "Counts"
  ) +
  scale_fill_manual(
    values = c("T" = "blue", "C" = "orange", "E" = "brown"),
    name = "Localization",labels = c("Total" = "Tissue", "C" = "cytosol", "E" = "ER")
  ) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

###{r Normalization for PCA and correlation plotting}
#CPM normalization
cpm <- apply(embr_SMG1i_df, 2, function(x) x/sum(as.numeric(x)) * 10^6)
#TMM normalization
norm.factors <- edgeR::calcNormFactors(cpm, method = "TMM")
norm <- edgeR::equalizeLibSizes(edgeR::DGEList(cpm, norm.factors = norm.factors))$pseudo.counts
#Spearman correlation test
correlationMatrix <- cor(normdata_ER, method = "spearman")
library(corrplot)
plot<-all2all(normdata_ER,cex=2)
#Plot
corrplot(correlationMatrix, order = 'hclust', 
         addrect = 2, addCoef.col = 'white', 
         number.cex = 0.7)
pairs(norm, pch = 19)

correlationMatrix <- cor(normdata_cyto, method = "spearman")
#Plot
corrplot(correlationMatrix, order = 'hclust', 
         addrect = 2, addCoef.col = 'white', 
         number.cex = 0.7,method="circle")
pairs(norm, pch = 19)
plot<-all2all(normdata_cyto,cex=2)

correlationMatrix <- cor(normdata, method = "spearman")
#Plot
corrplot(correlationMatrix, order = 'hclust', 
         addrect = 2, addCoef.col = 'white', 
         number.cex = 0.7)
pairs(norm, pch = 19)
plot<-all2all(normdata,cex=2)

############ PCA plot from riboSEq ER ########

ribodata <- read.table(file="/Volumes/AG_Landthaler/amanukyan/Projects/RiboSeq/Pipeline/run/sofia/output/FeatureCounts/cds_counts.txt", header = T, sep = "\t")
colnames(ribodata) <- lapply(colnames(ribodata), function(x) strsplit(x, split = "_S[0-9]+_")[[1]][1])
rownames(ribodata) <- ribodata[,1]
ribodata <- ribodata[,-1]
ribodata <- ribodata[,c(13:14,19,15,17,18)]
colnames(ribodata) <- gsub("output.bam_files.", "", colnames(ribodata))
colnames(ribodata) <- gsub(".bam", "", colnames(ribodata))


ribodataT <- read.table(file="/Volumes/AG_Landthaler/amanukyan/Projects/RiboSeq/Pipeline/run/sofia/output/FeatureCounts/cds_counts.txt", header = T, sep = "\t")
colnames(ribodataT) <- lapply(colnames(ribodataT), function(x) strsplit(x, split = "_S[0-9]+_")[[1]][1])
rownames(ribodataT) <- ribodataT[,1]
ribodataT <- ribodataT[,-1]
ribodataT <- ribodataT[,c(20:27)]
colnames(ribodataT) <- gsub("output.bam_files.", "", colnames(ribodataT))
colnames(ribodataT) <- gsub(".bam", "", colnames(ribodataT))


####

# remove low count genes
count_gene <- apply(ribodata, 1, max) 
ribodata <- ribodata[which(count_gene > 10),]

# control treatment
riboCond <- c("ctrl","ctrl","time 4 hours","time 4 hours","time 8 hours","time 8 hours")
riboCond<-factor(riboCond)

# ribo tpm
len_ribo <- len[match(rownames(ribodata), len$hgnc_symbol),]
len_ribo <- na.omit(len_ribo)
ribodata_tpm <- ribodata[intersect(rownames(ribodata), len_ribo$hgnc_symbol),]
ribodata_tpm <- sweep(ribodata_tpm * 1000, 1, len_ribo$len, FUN = "/")
ribodata_tpm <- sweep(ribodata_tpm * 1000000, 2, colSums(ribodata_tpm), FUN = "/")

####

# remove low count genes
count_geneT <- apply(ribodataT, 1, max) 
ribodataT <- ribodataT[which(count_geneT > 10),]

# control treatment
riboCondT <- c("ctrl","ctrl","ctrl","time 4 hours","time 4 hours","time 8 hours","time 8 hours","time 8 hours")
riboCondT<-factor(riboCondT)

# ribo tpm
len_riboT <- len[match(rownames(ribodataT), len$hgnc_symbol),]
len_riboT <- na.omit(len_riboT)
ribodata_tpmT <- ribodataT[intersect(rownames(ribodataT), len_riboT$hgnc_symbol),]
ribodata_tpmT <- sweep(ribodata_tpmT * 1000, 1, len_riboT$len, FUN = "/")
ribodata_tpmT <- sweep(ribodata_tpmT * 1000000, 2, colSums(ribodata_tpmT), FUN = "/")


library(debrowser)

normdata_RiboER <- getNormalizedMatrix(ribodata_tpm, method = "TMM")

normdata_pr_RiboER <- apply(normdata_RiboER,1,scale)
prtemp_RiboER <- prcomp(normdata_pr_RiboER)
normdata_pr_RiboER <- prtemp_RiboER$x
col<-colnames(normdata_pr_RiboER)
ggplot(data.frame(normdata_pr_RiboER), aes(x = PC1, y = PC2, color = riboCond)) + geom_point() + labs(title = "PCA ER RiboSeq") + geom_text(aes(label = col, vjust = 1.5, hjust = 0.5))


normdata_RiboT <- getNormalizedMatrix(ribodata_tpmT, method = "TMM")

normdata_pr_RiboT <- apply(normdata_RiboT,1,scale)
prtemp_RiboT <- prcomp(normdata_pr_RiboT)
normdata_pr_RiboT <- prtemp_RiboT$x
colT<-colnames(normdata_pr_RiboT)
ggplot(data.frame(normdata_pr_RiboT), aes(x = PC1, y = PC2, color = riboCondT)) + geom_point() + labs(title = "PCA Total RiboSeq") + geom_text(aes(label = colT, vjust = 1.5, hjust = 0.5))

correlationMatrixT <- cor(normdata_RiboT, method = "spearman")
correlationMatrixE <- cor(normdata_RiboER, method = "spearman")
correlationMatrixC <- cor(normdata_cyto, method = "spearman")
#Plot
corrplot(correlationMatrixT, order = 'hclust', 
         addrect = 2, addCoef.col = 'white', 
         number.cex = 0.7)
pairs(norm, pch = 19)

#Plot
corrplot(correlationMatrixE, order = 'hclust', 
         addrect = 2, addCoef.col = 'white', 
         number.cex = 0.7)
pairs(norm, pch = 19)
#Plot
corrplot(correlationMatrixC, order = 'hclust', 
         addrect = 2, addCoef.col = 'white', 
         number.cex = 0.7)
pairs(norm, pch = 19)


library(corrplot)
all2all(normdata_RiboER,cex=2)
all2all(normdata_RiboT,cex=2)
all2all(normdata_cyto,cex=2)
