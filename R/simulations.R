#' @export
simulations <- function(case_samples, control_samples, n = 100, methods = c("manova", "mlm", "mahdistmcd", "isoforest", "barbosa", "qn"), chr = NULL , start = NULL, end = NULL, epi_params = epi_parameters(),bump_cutoff = 0.1, min_cpg = 3, minoverlap =1000,  verbose = TRUE){
  
  if(is.null(case_samples)){
    stop("Please provide a valid 'case_samples'")
  }
  if(is.null(control_samples)){
    stop("Please provide a valid 'control_samples'")
  }
  total_control_samples <- ncol(control_samples)
  samples_ncol <- sample.int(total_control_samples, size = n, replace = FALSE)
  control_samples <- control_samples[,samples_ncol] 
  
  rst <- do.call(rbind, lapply(seq_len(length(methods)), function(i) {
    rst_case <- do.call(rbind, lapply(ncol(case_samples), function(ii) {
      epimutacions::epimutations(case_samples[,ii], 
                                 control_panel = control_samples, 
                                 method = methods[i], 
                                 chr = chr, 
                                 start = start, 
                                 end = end,  
                                 epi_params = epi_params,
                                 min_cpg = min_cpg,
                                 bump_cutoff = bump_cutoff,
                                 verbose = verbose)
    
      }))
  }))
  
  if(nrow(rst) == 0){
    rst <- NA
    return(rst)
  }
  gr <- GenomicRanges::makeGRangesFromDataFrame(as.data.frame(rst), seqnames.field = "chromosome", start.field = "start", end.field= "end")
  overlaps <- GenomicRanges::findOverlaps(gr,epi_validated, minoverlap = 0, maxgap = 600,  type ="equal")
  
  if(length(overlaps) != 0){
  keep_rst <- S4Vectors::queryHits(overlaps)
  keep_validated <- S4Vectors::subjectHits(overlaps)
  rst <- as.data.frame(rst[keep_rst,])
  validated <- as.data.frame(epi_validated[keep_validated,])
  results <- cbind(rst,validated[,1:4])
  results <- results[,c(1:7,12:15,8:11)]
  colnames(results) <- c("epi_id", "sample", "chromosome", "start", "end", "sz", "cpg_n", "chromosome_validated","start_validated","end_validated", "width_validated", "cpg_ids", "outlier_score", "outlier_significance", "outlier_direction")
  return(results)
  }else{
   return(NA)
  }
}
