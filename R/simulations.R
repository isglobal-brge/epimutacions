#' @export
simulations <- function(cohort, n = 100, methods = c("manova", "mlm", "mahdistmcd", "isoforest"), chr = NULL , start = NULL, end = NULL,  pvalue_cutoff = 0.01, outlier_score_cutoff = 0.6, min_cpg = 3, verbose = TRUE){
  
  if(is.null(cohort)){
    stop("Please provide a valid 'cohort'")
  }
  total_cohort_samples <- ncol(cohort)
  samples_ncol <- sample.int(total_cohort_samples, size = n, replace = FALSE)
  cohort_subset <- cohort[,samples_ncol] 
  
  rst <- do.call(rbind, lapply(seq_len(length(methods)), function(i) {
    rst <- do.call(rbind, lapply(seq_len(ncol(cohort_subset)), function(j){
      EpiMutations:: EpiMutations(case_sample = cohort_subset[,j], 
                                  control_panel = cohort_subset[,-j], 
                                  method = methods[i], 
                                  chr = chr, 
                                  start = start, 
                                  end = end,  
                                  pvalue_cutoff = pvalue_cutoff, 
                                  outlier_score_cutoff = outlier_score_cutoff,
                                  min_cpg = min_cpg,
                                  verbose = verbose)
    
      })) 
  }))
  if(nrow(rst) == 0){
    rst <- NA
   }
  return(rst)
}
