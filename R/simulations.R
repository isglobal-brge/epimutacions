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
    rst_case <- do.call(rbind, lapply(seq_len(ncol(case_samples)), function(ii) {
     epimutations(case_samples[,ii], 
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
  
  
  if(is.null(rst)){
    rst <- NA
    return(rst)
  }
  query <- GenomicRanges::makeGRangesFromDataFrame(as.data.frame(rst), seqnames.field = "chromosome", start.field = "start", end.field= "end")
  subject  <- epi_validated
  hits <- GenomicRanges::findOverlaps(query, subject, minoverlap = 0, maxgap = 600,  type ="equal")
  
  if(length(hits) != 0){
  overlaps <- GenomicRanges::pintersect(query[S4Vectors::queryHits(hits)], subject[S4Vectors::subjectHits(hits)])
  percentOverlap <- width(overlaps) / width(subject[S4Vectors::subjectHits(hits)])
  
  
  
  keep_rst <- S4Vectors::queryHits(hits)
  keep_validated <- S4Vectors::subjectHits(hits)
  rst <- as.data.frame(rst[keep_rst,])
  validated <- as.data.frame(epi_validated[keep_validated,])
  results <- cbind(rst,validated[,1:4])
  results$percent25 <- ifelse(percentOverlap >= 0.25, 1, 0)
  results$percent50 <- ifelse(percentOverlap >= 0.50, 1, 0)
  results$percent75 <- ifelse(percentOverlap >= 0.75, 1, 0)
  results$percent100 <- ifelse(percentOverlap >= 1, 1, 0)
  results$percent <- percentOverlap
  
  results <- results[,c(1:7,12:20)]
  colnames(results) <- c("epi_id", "sample", "chromosome", "start", "end", "sz", "cpg_n", "chromosome_validated","start_validated","end_validated", "width_validated", "percent25", "percent50", "percent75", "percent100", "percentpercent")
  return(results)
  }else{
    return(NA)
  }
  
 
}
