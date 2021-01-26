#' @title Epimutations analysis based on outlier detection methods
#' @description This function is similar to \link[epimutacions]{epimutations} with the
#' particularity that when is more than one case sample,  the remaining case samples
#' are included as controls. 
#' @param methy a GenomicRatioSet object containing the samples for the analysis. 
#' See the constructor function \link[minfi]{GenomicRatioSet}, \link[minfi]{makeGenomicRatioSetFromMatrix}. 
#' @param status a vector of 3 character string specifying the name of the status variable 
#' in the colData and the given name for controls and cases. These last 2 can be binomial. 
#' @param method a character string naming the outlier detection method to be used. 
#' This can be set as: \code{"manova"}, \code{"mlm"}, \code{"isoforest"}, \code{"mahdistmcd"}, 
#' \code{"barbosa"} and \code{"qn"}. 
#' The default is \code{"manova"}. 
#' For more information see \strong{Details}. 
#' @param chr a character string containing the sequence names to be analysed. The default value is \code{NULL}. 
#' @param start an integer specifying the start position. The default value is \code{NULL}.
#' @param end an integer specifying the end position. The default value is \code{NULL}.
#' @param epi_params The parameters for each method. See the function \link[epimutations]{epi_parameters}.  
#' @param bump_cutoff a numeric value of the estimate of the genomic profile above the 
#' cutoff or below the negative of the cutoff will be used as candidate regions. 
#' @param min_cpg an integer specifying the minimum CpGs number in a DMR.  
#' @param verbose logical. If TRUE additional details about the procedure will provide to the user. 
#' The default is TRUE. 
#' @details The function compares a case sample against a control panel to identify epimutations in the given 
#' sample. First, the DMRs are identified using the \link[bumphunter]{bumphunter} approach. 
#' After that, CpGs in those DMRs are tested in order to detect regions
#' with CpGs being outliers.  For that,  different anomaly detection methods can be selected:  
#'  * Multivariate Analysis of Variance (\code{"manova"}). \link[stats]{manova}
#'  * Multivariate Linear Model (\code{"mlm"})
#'  * Isolation Forest (\code{"isoforest"}) \link[isotree]{isolation.forest}
#'  * Robust Mahalanobis Distance (\code{"mahdistmcd"}) \link[robustbase]{covMcd}
#'  * Barbosa (\code{"barbosa"})
#'  * Qn (\code{"Qn"})
#' @return The function returns an object of class tibble containing the outliers regions.  
#' The results are composed by the following columns: 
#' * \code{epi_id}: the name of the anomaly detection method that has been used to detect the epimutation
#' * \code{sample}: the name of the sample where the epimutation was found.
#' * \code{chromosome}, \code{start} and \code{end}: the region coordinates.
#' * \code{sz}: the number of base pairs in the region.
#' * \code{cpg_n}: number of CpGs in the region.
#' * \code{cpg_ids}: differentially methylated CpGs names.
#' * \code{outlier_score}: the outlier score of the DMRs. 
#' The outlier score only is available in the output of manova, mlm, isolation forest and qn. 
#' * \code{outlier_significance}: the outlier significance of the region.  Only it is 
#'  available for manova and mlm. 
#' * \code{outlier_direction}: describes if the epimutation is hypermethylated or hypomethylated. 
#' The outlier direction is given by the methods barbosa, qn, isoforest, manova and mlm.
#' @examples 
#' \dontrun{
#' library(epimutacions)
#' data(methy)
#' 
#' #Find epimutations in GSM2562701 sample of methy dataset
#' 
#' epimutations_one_leave_out(methy)
#' }
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

