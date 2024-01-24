## load library
library(ribosomeProfilingQC)
library(AnnotationDbi)
library(Rsamtools)
library(BSgenome.Hsapiens.UCSC.hg38)
#prepare cds annottaions
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(Rsamtools)
library(gridExtra)
library(GenomicFeatures)
library(reshape2)
library(ggplot2)

# genome and taxonomy database
genome <- Hsapiens
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene ## give it a short name
CDS <- prepareCDS(txdb)

setwd("/Volumes/AG_Landthaler/sofia/TOTAL")

# Path to the directory containing BAM files
bam_files <- list.files(pattern = "_er.bam$|_cytosol.bam$", full.names = FALSE)

# make coverage and metaplot
cvgs.utr5 <- coverageDepth(RPFs = bam_files, gtf = txdb, region="utr5")
cvgs.CDS <- coverageDepth(RPFs = bam_files, gtf = txdb, region="cds")
cvgs.utr3 <- coverageDepth(RPFs = bam_files, gtf = txdb, region="utr3")
save(cvgs.utr5, cvgs.CDS, cvgs.utr3, file = "cvgs.RData")
load("cvgs.RData")

# cytosol and ER of one sample
g <- list()
for(i in seq(1, length(bam_files), 2)){
  title_bam <- strsplit(bam_files[i], split = "\\.bam")[[1]][1]
  print(title_bam)
  plotdata1 <- metaPlot(cvgs.utr5, cvgs.CDS, cvgs.utr3, plot=TRUE, sample = i)
  plotdata2 <- metaPlot(cvgs.utr5, cvgs.CDS, cvgs.utr3, plot=TRUE, sample = i+1)
  plotdata1_v <- unlist(plotdata1, use.names = TRUE)
  plotdata2_v <- unlist(plotdata2, use.names = TRUE)
  plotdata_merged <- rbind(plotdata1_v, plotdata2_v)
  plotdata_merged <- plotdata_merged[,1:200]
  plotdata_merged <- t(plotdata_merged)
  colnames(plotdata_merged) <- c("cytosol", "ER")
  plotdata_merged <- melt(plotdata_merged)
  colnames(plotdata_merged) <- c("Region", "Location", "value")
  g[[i]] <- ggplot2::ggplot(plotdata_merged, aes_string(x = "Region", y = "value", color = "Location", group = "Location")) + 
    geom_line() + scale_x_discrete(breaks = plotdata_merged$Region[100], labels = "START") + geom_vline(xintercept = 100) + theme_classic() + 
    labs(title = title_bam) + ylim(0,0.2)
}
g2 <- g[!unlist(lapply(g,is.null))]
ggpubr::ggarrange(plotlist = g2, ncol = 5, nrow = 5, common.legend = TRUE)