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

setwd("/Volumes/AG_Landthaler/sofia/cytosol_membrane_bam")
source("Artur_function_coverage_100.R")
# Path to the directory containing BAM files
bam_files <- list.files(pattern = "_er.bam$|_cytosol.bam$", full.names = FALSE)
# make coverage and metaplot
cvgs.utr5 <- coverageDepth(RPFs = bam_files, gtf = txdb, region="utr5")
cvgs.CDS <- coverageDepth(RPFs = bam_files, gtf = txdb, region="cds")
cvgs.utr3 <- coverageDepth(RPFs = bam_files, gtf = txdb, region="utr3")
save(cvgs.utr5, cvgs.CDS, cvgs.utr3, file = "cvgs.RData")
load("cvgs.RData")

#bam_files_er[3]

# ER samples
bam_files_er <- bam_files[grepl("^E[0-9]R2", bam_files)]
bam_files_er <- bam_files_er[grepl("er.bam$", bam_files_er)]
ind1 <- which(bam_files == bam_files_er[1])
ind2 <- which(bam_files == bam_files_er[2])
ind3 <- which(bam_files == bam_files_er[3])
plotdata1 <- metaPlotNT(cvgs.utr5, cvgs.CDS, sample = ind1)
plotdata2 <- metaPlotNT(cvgs.utr5, cvgs.CDS, sample = ind2)
plotdata3 <- metaPlotNT(cvgs.utr5, cvgs.CDS, sample = ind3)
plotdata_merged <- rbind(plotdata1, plotdata2, plotdata3)
plotdata_merged <- t(plotdata_merged)
colnames(plotdata_merged) <- c("Ctrl", "T4", "T8")
plotdata_merged <- melt(plotdata_merged)
colnames(plotdata_merged) <- c("Condition", "Location", "value")
labels <- c("START")
ggplot2::ggplot(plotdata_merged, aes_string(x = "Condition", y = "value", color = "Location", group = "Location")) + 
  geom_line() + scale_x_discrete(breaks = plotdata_merged$Condition[34], labels = labels) + geom_vline(xintercept = 33) + theme_classic()

# Total samples
bam_files_total <- bam_files[grepl("^TT[0-9]R2", bam_files)]
bam_files_total <- bam_files_total[grepl("er.bam$", bam_files_total)]
ind1 <- which(bam_files == bam_files_total[1])
ind2 <- which(bam_files == bam_files_total[2])
ind3 <- which(bam_files == bam_files_total[3])
plotdata1 <- metaPlotNT(cvgs.utr5, cvgs.CDS, sample = ind1)
plotdata2 <- metaPlotNT(cvgs.utr5, cvgs.CDS, sample = ind2)
plotdata3 <- metaPlotNT(cvgs.utr5, cvgs.CDS, sample = ind3)
plotdata_merged <- rbind(plotdata1, plotdata2, plotdata3)
plotdata_merged <- t(plotdata_merged)
colnames(plotdata_merged) <- c("Ctrl", "T4", "T8")
plotdata_merged <- melt(plotdata_merged)
colnames(plotdata_merged) <- c("Condition", "Location", "value")
labels <- c("START")
ggplot2::ggplot(plotdata_merged, aes_string(x = "Condition", y = "value", color = "Location", group = "Location")) + 
  geom_line() + scale_x_discrete(breaks = plotdata_merged$Condition[34], labels = labels) + geom_vline(xintercept = 33) + theme_classic()


#SP samples###

setwd("/Volumes/AG_Landthaler/sofia/sp_bam")
# Path to the directory containing BAM files
bam_files <- list.files(path = "/Volumes/AG_Landthaler/sofia/sp_bam",pattern = "_sp.bam$", full.names = FALSE)
# make coverage and metaplot
cvgs.utr5 <- coverageDepth(RPFs = bam_files, gtf = txdb, region="utr5")
cvgs.CDS <- coverageDepth(RPFs = bam_files, gtf = txdb, region="cds")
cvgs.utr3 <- coverageDepth(RPFs = bam_files, gtf = txdb, region="utr3")
save(cvgs.utr5, cvgs.CDS, cvgs.utr3, file = "cvgs_sp.RData")
load("cvgs_sp.RData")
# ER samples
bam_files_er <- bam_files[grepl("^E[0-9]R2", bam_files)]
bam_files_er <- bam_files_er[grepl("sp.bam$", bam_files_er)]
ind1 <- which(bam_files == bam_files_er[1])
ind2 <- which(bam_files == bam_files_er[2])
ind3 <- which(bam_files == bam_files_er[3])
plotdata1 <- metaPlotNT(cvgs.utr5, cvgs.CDS, sample = ind1)
plotdata2 <- metaPlotNT(cvgs.utr5, cvgs.CDS, sample = ind2)
plotdata3 <- metaPlotNT(cvgs.utr5, cvgs.CDS, sample = ind3)
plotdata_merged <- rbind(plotdata1, plotdata2, plotdata3)
plotdata_merged <- t(plotdata_merged)
colnames(plotdata_merged) <- c("Ctrl", "T4", "T8")
plotdata_merged <- melt(plotdata_merged)
colnames(plotdata_merged) <- c("Condition", "Location", "value")
labels <- c("START")
ggplot2::ggplot(plotdata_merged, aes_string(x = "Condition", y = "value", color = "Location", group = "Location")) + 
  geom_line() + scale_x_discrete(breaks = plotdata_merged$Condition[34], labels = labels) + geom_vline(xintercept = 33) + theme_classic()
# Total samples
bam_files_total <- bam_files[grepl("^TT[0-9]R2", bam_files)]
bam_files_er <- bam_files_er[grepl("sp.bam$", bam_files_er)]
ind1 <- which(bam_files == bam_files_total[1])
ind2 <- which(bam_files == bam_files_total[2])
ind3 <- which(bam_files == bam_files_total[3])
plotdata1 <- metaPlotNT(cvgs.utr5, cvgs.CDS, sample = ind1)
plotdata2 <- metaPlotNT(cvgs.utr5, cvgs.CDS, sample = ind2)
plotdata3 <- metaPlotNT(cvgs.utr5, cvgs.CDS, sample = ind3)
plotdata_merged <- rbind(plotdata1, plotdata2, plotdata3)
plotdata_merged <- t(plotdata_merged)
colnames(plotdata_merged) <- c("Ctrl", "T4", "T8")
plotdata_merged <- melt(plotdata_merged)
colnames(plotdata_merged) <- c("Condition", "Location", "value")
labels <- c("START")
ggplot2::ggplot(plotdata_merged, aes_string(x = "Condition", y = "value", color = "Location", group = "Location")) + 
  geom_line() + scale_x_discrete(breaks = plotdata_merged$Condition[34], labels = labels) + geom_vline(xintercept = 33) + theme_classic()


plotTranscriptc <- function(reads, tx_name,
                           col=c("Frame_0" = "#009E73",
                                 "Frame_1" = "#D55E00",
                                 "Frame_2" = "#0072B2")){
  stopifnot(is(reads, "GRanges"))
  if(length(reads$tx_name)!=length(reads) ||
     length(reads$position)!=length(reads) ||
     length(reads$posToStop)!=length(reads) ||
     length(reads$readingFrame)!=length(reads) ||
     length(reads$gene_id)!=length(reads)){
    stop("reads must be a result of assignReadingFrame")
  }
  l <- length(tx_name)
  op <- par(mfrow = c(ceiling(l/floor(sqrt(l))), floor(sqrt(l))))
  on.exit(par(op))
  heights <- list()
  for(i in tx_name){
    x.sub <- reads[reads$tx_name %in% i]
    if(length(x.sub)<1){
      adist <- adist(reads$tx_name, i)
      adist.min <- min(adist, na.rm = TRUE)
      possibleTxName <- 
        unique(reads$tx_name[adist==adist.min & 
                               !is.na(reads$tx_name)])
      warning("No reads in ", i, ".",
              "The closest transcripts are: ",
              paste(possibleTxName, collapse = ", "))
      heights[[i]] <- numeric()
    }else{
      d <- table(mcols(x.sub)[, c("readingFrame", "position")])
      CDS.size <- x.sub[1]$position + x.sub[1]$posToStop + 3
      at <- as.character(seq.int(CDS.size))
      height <- colSums(d)
      height <- height[at]
      height[is.na(height)] <- 0
      names(height) <- at
      cols <- col[apply(d, 2, function(.ele) which(.ele!=0))]
      names(cols) <- colnames(d)
      cols <- cols[at]
      names(cols) <- at
      barplot(xlim=c(0, CDS.size),ylim=c(0,100),# change ylim
              height=height,
              xlab = paste("Base position relative to", i,"CDS"),
              ylab = "Number of reads",
              col = cols, border=cols,  main = i)
      legend("topleft", legend = names(col),
             fill = col, border = col, bg = NA, box.col = NA) 
      heights[[i]] <- height
    }
  }
  return(invisible(heights))
}

#####ER

#t0
bamfilename_0="/Volumes/AG_Landthaler/amanukyan/Projects/RiboSeq/Pipeline/run/sofia/output/bam_files/E0R2_RPI10_S9.bam"
yieldSize <- 10000000
bamfile_0 <- BamFile(bamfilename_0, yieldSize = yieldSize)
pc_0 <- getPsiteCoordinates(bamfile_0, bestpsite = 13)
pc.sub_0 <- pc_0[pc_0$qwidth %in% c(27, 28)]
pc.sub_0 <- assignReadingFrame(pc.sub_0, CDS)
plotTranscriptc(pc.sub_0, c("ENST00000221566.7","ENST00000675538.1","ENST00000324460.7","ENST00000620073.4","ENST00000369085.8","ENST00000227507.3"))
#t4
bamfilename_4="/Volumes/AG_Landthaler/amanukyan/Projects/RiboSeq/Pipeline/run/sofia/output/bam_files/E4R2_RPI12_S11.bam"
yieldSize <- 10000000
bamfile_4 <- BamFile(bamfilename_4, yieldSize = yieldSize)
pc_4 <- getPsiteCoordinates(bamfile_4, bestpsite = 13)
pc.sub_4 <- pc_4[pc_4$qwidth %in% c(27, 28)]
pc.sub_4 <- assignReadingFrame(pc.sub_4, CDS)
plotTranscriptc(pc.sub_4, c("ENST00000221566.7","ENST00000675538.1","ENST00000324460.7","ENST00000620073.4","ENST00000369085.8","ENST00000227507.3")) 
#t8
bamfilename_8="/Volumes/AG_Landthaler/amanukyan/Projects/RiboSeq/Pipeline/run/sofia/output/bam_files/E8R2_RPI14_S13.bam"
yieldSize <- 10000000
bamfile_8 <- BamFile(bamfilename_8, yieldSize = yieldSize)
pc_8 <- getPsiteCoordinates(bamfile_8, bestpsite = 13)
pc.sub_8 <- pc_8[pc_8$qwidth %in% c(27, 28)]
pc.sub_8 <- assignReadingFrame(pc.sub_8, CDS)
plotTranscriptc(pc.sub_8, c("ENST00000677562.1","ENST00000221566.7","ENST00000675538.1","ENST00000324460.7","ENST00000620073.4","ENST00000369085.8","ENST00000227507.3")) 

plotTranscriptc(pc.sub_0, c("ENST00000221566.7","ENST00000679138.1"))
plotTranscriptc(pc.sub_4, c("ENST00000221566.7","ENST00000677562.1"))
plotTranscriptc(pc.sub_8, c("ENST00000221566.7","ENST00000677562.1"))
####Total

#t0
bamfilename_0t="/Volumes/AG_Landthaler/amanukyan/Projects/RiboSeq/Pipeline/run/sofia/output/bam_files/TT0R2_RPI15_S14.bam"
yieldSize <- 10000000
bamfile_0t <- BamFile(bamfilename_0t, yieldSize = yieldSize)
pc_0t <- getPsiteCoordinates(bamfile_0t, bestpsite = 13)
pc.sub_0t <- pc_0t[pc_0t$qwidth %in% c(27, 28)]
pc.sub_0t <- assignReadingFrame(pc.sub_0t, CDS)
plotTranscriptc(pc.sub_0t, c("ENST00000675538.1","ENST00000324460.7","ENST00000620073.4","ENST00000369085.8","ENST00000680688.1"))
#t4
bamfilename_4t="/Volumes/AG_Landthaler/amanukyan/Projects/RiboSeq/Pipeline/run/sofia/output/bam_files/TT4R2_RPI16_S15.bam"
bamfile_4t <- BamFile(bamfilename_4t, yieldSize = yieldSize)
pc_4t <- getPsiteCoordinates(bamfile_4t, bestpsite = 13)
pc.sub_4t<- pc_4t[pc_4t$qwidth %in% c(27, 28)]
pc.sub_4t <- assignReadingFrame(pc.sub_4t, CDS)
plotTranscriptc(pc.sub_4t, c("ENST00000675538.1","ENST00000324460.7","ENST00000620073.4","ENST00000369085.8","ENST00000680688.1"))
#t8
bamfilename_8t="/Volumes/AG_Landthaler/amanukyan/Projects/RiboSeq/Pipeline/run/sofia/output/bam_files/TT8R2_RPI17_S16.bam"
bamfile_8t<- BamFile(bamfilename_8t, yieldSize = yieldSize)
pc_8t <- getPsiteCoordinates(bamfile_8t, bestpsite = 13)
pc.sub_8t <- pc_8t[pc_8t$qwidth %in% c(27, 28)]
pc.sub_8t <- assignReadingFrame(pc.sub_8t, CDS)
plotTranscriptc(pc.sub_8t, c("ENST00000675538.1","ENST00000324460.7","ENST00000620073.4","ENST00000369085.8","ENST00000680688.1"))

