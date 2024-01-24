metaPlotNT <- function (UTR5coverage, CDScoverage, sample, cover = c(UTR5 = 100, CDS = 300), ...) 
{
  xaxis <- "RPFs"
  if (!is.list(UTR5coverage) || !is.list(CDScoverage)) {
    stop("UTR5coverage, CDScoverage and UTR3coverage must be\n         output of coverageDepth")
  }
  if (!xaxis %in% names(UTR5coverage) || !xaxis %in% names(CDScoverage)) {
    stop("UTR5coverage, CDScoverage and UTR3coverage must be \n output of coverageDepth.", "And RPFs or mRNA must be available.")
  }
  if (length(sample) > 1) {
    stop("Only one sample are computed")
  }
  if (!all(c("UTR5", "CDS") %in% names(cover))) {
    stop("Covers should be a integer vector with names UTR5, CDS and UTR3.")
  }
  regions <- list(UTR5 = UTR5coverage[[xaxis]][["coverage"]][[sample]], 
                  CDS = CDScoverage[[xaxis]][["coverage"]][[sample]])
  cover <- cover[names(regions)]
  regions <- mapply(regions, cover, FUN = function(.ele, .len) {
    .ele[lengths(.ele) >= .len]
  }, SIMPLIFY = FALSE)
  ids <- Reduce(intersect, lapply(regions, names))
  regions <- lapply(regions, `[`, i = ids)
  metagene <- mapply(regions, cover, FUN = function(.ele, .len) {
    if(.len == cover["CDS"]) {
      .ele <- sapply(.ele, function(x) x[1:.len])
    } else {
      .ele <- sapply(.ele, function(x) tail(x,.len))
    }
    .ele <- RleList(.ele)
    l <- rep(.len, length(.ele))
    
    # ir <- IRanges(start = 1:.len, end = 1:.len, width = rep(1,.len))
    # ir <- IRangesList(rep(list(ir), length(.ele)))
    # vws <- Views(.ele, ir)
    
    ir <- IRanges(1, width = l)
    print(ir)
    # tile <- tile(ir, n = .len/5)
    tile <- tile(ir, width = 3)
    names(tile) <- names(.ele)
    vws <- Views(.ele, tile)
    
    vms <- viewMeans(vws)
    vms <- do.call(rbind, as.list(vms))
    vms <- colMeans(vms)
    vms
  }, SIMPLIFY = FALSE)
  metagene <- unlist(metagene, use.names = TRUE)
  return(invisible(metagene))
}
