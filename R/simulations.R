#' @export
simulations <- function(cohort, status, fd, n = 100, methods = c("manova", "mlm", "mahdistmcd", "isoforest", "barbosa", "qn"), chr = NULL , start = NULL, end = NULL, epi_params = epi_parameters(),bump_cutoff = 0.1, min_cpg = 3, minoverlap =1000,  verbose = TRUE){
  
  if(is.null(cohort)){
    stop("Please provide a valid 'cohort'")
  }
  total_cohort_samples <- ncol(cohort)
  samples_ncol <- sample.int(total_cohort_samples, size = n, replace = FALSE)
  cohort <- cohort[,samples_ncol] 
  status <- status[samples_ncol]
  control_panel <- cohort[,status == 0]
  case_samples <- cohort[,status == 1]
  
  rst <- do.call(rbind, lapply(seq_len(length(methods)), function(i) {
      EpiMutations::EpiMutations(case_samples = case_samples, 
                                  control_panel = control_panel, 
                                  fd = fd, 
                                  method = methods[i], 
                                  chr = chr, 
                                  start = 1, 
                                  end = 1000000,  
                                  epi_params = epi_params,
                                  min_cpg = min_cpg,
                                  bump_cutoff = bump_cutoff,
                                  verbose = verbose)
    
      })) 
  
  if(nrow(rst) == 0){
    rst <- NA
    return(rst)
  }
  gr <- GenomicRanges::makeGRangesFromDataFrame(as.data.frame(rst), seqnames.field = "chromosome", start.field = "start", end.field= "end")
  keep <- GenomicRanges::findOverlaps(gr,epi_validated, minoverlap = 0, maxgap= 100000000,  type ="equal")
  keep <- S4Vectors::queryHits(keep)
  rst <- rst[keep,]
  
  methods <- strsplit(rst$epi_id, split = "_")
  methods <- lapply(methods,function(x)  x[2])
  methods <- unlist(methods)
  rst$epi_id <- methods
  method <- unique(methods)
  
  total <- NULL
  chr5 <- NULL
  chr17 <- NULL
  chr19 <- NULL
  
  for(i in 1:length(method)){
    total <- c(total, length(which(rst$epi_id ==method[i]))/10000)
    chr5 <- c(chr5, nrow(rst[which(rst$epi_id ==method[i] & rst$chromosome == "chr5"),])/10000)
    chr17 <- c(chr17,nrow(rst[which(rst$epi_id ==method[i] & rst$chromosome == "chr17"),])/10000)
    chr19 <- c(chr19,nrow(rst[which(rst$epi_id ==method[i] & rst$chromosome == "chr19"),])/10000) 
  }
  
  results <- as.data.frame(cbind(total,chr5,chr17,chr19))
  rownames(results) <- method
  
  return(results)
}
