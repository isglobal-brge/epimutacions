#' @export 
epimutations_one_leave_out <- function(methy, status = c("status", "control", "case"), method = "manova", chr = NULL, start = NULL, end = NULL, epi_params = epi_parameters(), bump_cutoff =  0.1, min_cpg = 3, verbose = TRUE){
  
  if(class(methy) != "GenomicRatioSet")
  {
    stop("Input data 'methy' must be a 'GenomicRatioSet'. 
         'makeGenomicRatioSetFromMatrix' function from 'minfi' package 
         can be useful to create a 'GenomicRatioSet' class object")
  }
  
  pd <- as.data.frame(SummarizedExperiment::colData(methy))
  
  if(length(status) != 3){
    stop("'status' must specify (1) the colData column, (2) the control level and (3) cases level names")
  }
  if(!status[1] %in% colnames(pd)){
    stop("The variable name '", status[1], "' is not in 'colData'")
  }
  if(!status[2] %in% unique(pd[,status[1]])){
    stop(" '", status[2], "' is not a level of '", status[1], "'")
    
  }
  if(!status[3] %in% unique(pd[,status[1]])){
    stop(" '", status[3], "' is not a level of '", status[1], "'")
    
  }
  
  case_samples <- pd[, status[1]] == status[3]
  case_samples_names <- colnames(methy[,case_samples])
  rm(case_samples)
  
  rst <- do.call(rbind, lapply(case_samples_names, function(case) {
    case_samples <- methy[,case]
    control_panel <- methy[,-which(colnames(methy) == case)]
    epimutacions::epimutations(case_samples, control_panel, method, chr, start, end, epi_params, bump_cutoff, min_cpg, verbose)
  }))
  return(rst)
}

