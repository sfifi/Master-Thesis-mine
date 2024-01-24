## load library
library(ribosomeProfilingQC)
library(AnnotationDbi)
library(Rsamtools)
#load genome
library(BSgenome.Hsapiens.UCSC.hg38)
#prepare cds annottaions
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(Rsamtools)
library(gridExtra)
library(GenomicFeatures)

genome <- Hsapiens
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene ## give it a short name
CDS <- prepareCDS(txdb)

processBAM <- function(bamfilename, txdb, CDS, genome) {
  # Index BAM
  indexBam(bamfilename)
  # Open BAM file
  yieldSize <- 10000000
  bamfile <- BamFile(bamfilename, yieldSize = yieldSize)
  # Estimate P-site
  estimatePsite(bamfile, CDS, genome)
  estimatePsite(bamfile, CDS, genome, anchor = "3end")
  # Calculate coverage depth for different regions
  cvgs.utr5 <- coverageDepth(RPFs = bamfilename, gtf = txdb, region="utr5")
  cvgs.CDS <- coverageDepth(RPFs = bamfilename, gtf = txdb, region="cds")
  cvgs.utr3 <- coverageDepth(RPFs = bamfilename, gtf = txdb, region="utr3")
  # Determine ylim based on filename keywords
  if (grepl("cytosol", bamfilename, ignore.case = TRUE)) {
    y_limit <- c(0, 0.35)
  } else if (grepl("er", bamfilename, ignore.case = TRUE)) {
    y_limit <- c(0, 0.15)
  } else {
    y_limit <- c(0, 0.4) # Default ylim as in your original code
  }
  
  # Generate and return meta plot
  plot <- metaPlot(cvgs.utr5, cvgs.CDS, cvgs.utr3, sample=1, plot=FALSE, ylim = y_limit)
  # Extract the desired portion of the filename for the title
  plot_title <- sub("_(.*)", "", basename(bamfilename))
  
  # Add "er" or "cytosol" to the title if they appear in the BAM filename
  if (grepl("cytosol", bamfilename, ignore.case = TRUE)) {
    plot_title <- paste(plot_title, "cytosol")
  } else if (grepl("er", bamfilename, ignore.case = TRUE)) {
    plot_title <- paste(plot_title, "er")
  }
  # Add the constructed title to the plot
  title(main = plot_title)
  cat(paste("Processed:", bamfilename, "\n"))
  return(plot)
}
# Path to the directory containing BAM files
bam_directory <- "/Volumes/AG_Landthaler/sofia/ER/"
# List all BAM files in the directory
bam_files <- list.files(path = bam_directory, pattern = "\\.bam$", full.names = TRUE)
# Apply the function to each BAM file and collect plots
plots <- lapply(bam_files, processBAM, txdb=txdb, CDS=CDS, genome=genome)
#single plots
single_plot <- processBAM(bam_files["x"], txdb=txdb, CDS=CDS, genome=genome)
print(single_plot)


############### SP bam genes ###############
processBAM <- function(bamfilename, txdb, CDS, genome) {
  # Index BAM
  indexBam(bamfilename)
  # Open BAM file
  yieldSize <- 10000000
  bamfile <- BamFile(bamfilename, yieldSize = yieldSize)
  # Estimate P-site
  estimatePsite(bamfile, CDS, genome)
  estimatePsite(bamfile, CDS, genome, anchor = "3end")
  # Calculate coverage depth for different regions
  cvgs.utr5 <- coverageDepth(RPFs = bamfilename, gtf = txdb, region="utr5")
  cvgs.CDS <- coverageDepth(RPFs = bamfilename, gtf = txdb, region="cds")
  cvgs.utr3 <- coverageDepth(RPFs = bamfilename, gtf = txdb, region="utr3")
  # Determine ylim based on filename keywords
  if (grepl("cytosol", bamfilename, ignore.case = TRUE)) {
    y_limit <- c(0, 0.35)
  } else if (grepl("er", bamfilename, ignore.case = TRUE)) {
    y_limit <- c(0, 0.15)
  } else {
    y_limit <- c(0, 0.4) # Default ylim as in your original code
  }
  
  # Generate and return meta plot
  plot <- metaPlot(cvgs.utr5, cvgs.CDS, cvgs.utr3, sample=1, plot=FALSE, ylim = y_limit)
  # Extract the desired portion of the filename for the title
  plot_title <- sub("_(.*)", "", basename(bamfilename))
  
  # Add "er" or "cytosol" to the title if they appear in the BAM filename
  if (grepl("cytosol", bamfilename, ignore.case = TRUE)) {
    plot_title <- paste(plot_title, "cytosol")
  } else if (grepl("er", bamfilename, ignore.case = TRUE)) {
    plot_title <- paste(plot_title, "er")
  }
  # Add the constructed title to the plot
  title(main = plot_title)
  cat(paste("Processed:", bamfilename, "\n"))
  return(plot)
}
# Path to the directory containing BAM files
bam_directory <- "/Volumes/AG_Landthaler/sofia/ER/"
# List all BAM files in the directory
bam_files <- list.files(path = bam_directory, pattern = "\\.bam$", full.names = TRUE)
# Apply the function to each BAM file and collect plots
plots <- lapply(bam_files, processBAM, txdb=txdb, CDS=CDS, genome=genome)
#single plots
single_plot <- processBAM(bam_files["x"], txdb=txdb, CDS=CDS, genome=genome)
print(single_plot)


processBAM_SP <- function(bamfilename, txdb, CDS, genome) {
  # Index BAM
  indexBam(bamfilename)
  # Open BAM file
  yieldSize <- 10000000
  bamfile <- BamFile(bamfilename, yieldSize = yieldSize)
  # Estimate P-site
  estimatePsite(bamfile, CDS, genome)
  estimatePsite(bamfile, CDS, genome, anchor = "3end")
  # Calculate coverage depth for different regions
  cvgs.utr5 <- coverageDepth(RPFs = bamfilename, gtf = txdb, region="utr5")
  cvgs.CDS <- coverageDepth(RPFs = bamfilename, gtf = txdb, region="cds")
  cvgs.utr3 <- coverageDepth(RPFs = bamfilename, gtf = txdb, region="utr3")
  
  # Generate meta plot without specifying ylim
  plot <- metaPlot(cvgs.utr5, cvgs.CDS, cvgs.utr3, sample=1, plot=FALSE,ylim(0,0.1))
  
  # Extract the desired portion of the filename for the title
  plot_title <- sub("_(.*)", "", basename(bamfilename))
  # Add the constructed title to the plot
  title(main = plot_title)
  
  # Set color based on sample name
  if (grepl("E", bamfilename, ignore.case = TRUE)) {
    col <- "red"
  } else if (grepl("C", bamfilename, ignore.case = TRUE)) {
    col <- "blue"
  } else if (grepl("T", bamfilename, ignore.case = TRUE)) {
    col <- "green"
  } else {
    col <- "black" # Default color
  }
  plot$col <- col
  
  cat(paste("Processed:", bamfilename, "\n"))
  return(plot)
}

# Path to the directory containing BAM files
bam_directory <- "/Volumes/AG_Landthaler/sofia/sp_bam"
# List all BAM files in the directory
bam_files <- list.files(path = bam_directory, pattern = "\\.bam$", full.names = TRUE)
# Apply the function to each BAM file and collect plots
plots <- lapply(bam_files, processBAM_SP, txdb=txdb, CDS=CDS, genome=genome)
#Example single plot
#single_plot <- processBAM_SP(bam_files[10], txdb=txdb, CDS=CDS, genome=genome)
#print(single_plot)

for (p in plots) {
  print(p)
}

##### try to see individual genes ######

bamfilename="/Volumes/AG_Landthaler/sofia/ER/E0R2_RPI10_S9.bam.assigned_er.bam"
#estimatePsite(bamfile, CDS, genome) #it is important to add in bestpsite what is obtained here
pc <- getPsiteCoordinates(bamfile, bestpsite = 13)
## for this QC demo, we will only use reads length of 28-29 nt.
pc.sub <- pc[pc$qwidth %in% c(27,28)]
strandPlot(pc.sub, CDS)
pc.sub <- readsDistribution(pc.sub, txdb, las=2)
pc.sub <- assignReadingFrame(pc.sub, CDS)
plotDistance2Codon(pc.sub)
plotTranscript(pc.sub, c("ENST00000674235.1"))
cvg <- frameCounts(pc.sub, coverageRate=TRUE,level="tx")
ORFscore <- getORFscore(pc.sub)
## Following code will plot the ORFscores vs coverage.
## Try it by removing the '#'. 
plot(cvg[names(ORFscore)], ORFscore,
    xlab="coverage ORF1", ylab="ORF score",
    type="p", pch=16, cex=.5, xlim=c(0, 1))

bamfilename="/Volumes/AG_Landthaler/sofia/ER/E4R2_RPI12_S11.bam.assigned_er.bam"
yieldSize <- 10000000
bamfile <- BamFile(bamfilename, yieldSize = yieldSize)
estimatePsite(bamfile, CDS, genome)
pc <- getPsiteCoordinates(bamfile, bestpsite = 13)

pc <- getPsiteCoordinates(bamfile, 13)
pc.sub <- pc[pc$qwidth %in% c(27, 28)]
pc.sub <- assignReadingFrame(pc.sub, CDS)
plotDistance2Codon(pc.sub)

## for this QC demo, we will only use reads length of 28-29 nt.
pc.sub <- pc[pc$qwidth %in% c(27,28)]
strandPlot(pc.sub, CDS)
pc.sub <- readsDistribution(pc.sub, txdb, las=2)
pc.sub <- assignReadingFrame(pc.sub, CDS)
plotDistance2Codon(pc.sub)
plotTranscript(pc.sub, c("ENST00000324460.7"))
bamfilename="/Volumes/AG_Landthaler/sofia/ER//E8R2_RPI14_S13.bam.assigned_er.bam"
yieldSize <- 10000000
bamfile <- BamFile(bamfilename, yieldSize = yieldSize)
estimatePsite(bamfile, CDS, genome)
pc <- getPsiteCoordinates(bamfile, bestpsite = 13)
## for this QC demo, we will only use reads length of 28-29 nt.
pc.sub <- pc[pc$qwidth %in% c(27,28)]
strandPlot(pc.sub, CDS)
pc.sub <- readsDistribution(pc.sub, txdb, las=2)
pc.sub <- assignReadingFrame(pc.sub, CDS)
plotDistance2Codon(pc.sub)
plotTranscript(pc.sub, c("ENST00000672860.3"))

bamfilename="/Volumes/AG_Landthaler/sofia/ER/E0R2_RPI10_S9.bam.assigned_er.bam"
yieldSize <- 10000000
bamfile <- BamFile(bamfilename, yieldSize = yieldSize)
estimatePsite(bamfile, CDS, genome)
pc <- getPsiteCoordinates(bamfile, bestpsite = 13)
## for this QC demo, we will only use reads length of 28-29 nt.
pc.sub <- pc[pc$qwidth %in% c(27,28)]
strandPlot(pc.sub, CDS)
pc.sub <- readsDistribution(pc.sub, txdb, las=2,precedence="CDS")
pc.sub <- assignReadingFrame(pc.sub, CDS)
plotDistance2Codon(pc.sub)
plotTranscript(pc.sub, c("ENST00000344347.6"))







