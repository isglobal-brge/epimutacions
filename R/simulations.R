#' @export
simulations <- function(cohort, status, fd, n = 100, methods = c("manova", "mlm", "mahdistmcd", "isoforest", "barbosa", "qn"), chr = NULL , start = NULL, end = NULL, epi_params = epi_params(),bump_cutoff = 0.1, min_cpg = 3, verbose = TRUE){
  
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
                                  start = start, 
                                  end = end,  
                                  epi_params = epi_params,
                                  min_cpg = min_cpg,
                                  bump_cutoff = bump_cutoff,
                                  verbose = verbose)
    
      })) 
  if(nrow(rst) == 0){
    rst <- NA
   }
  return(rst)
}
