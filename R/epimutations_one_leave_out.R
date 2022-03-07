#' @title Epimutations analysis based on outlier detection methods
#' @description This function is similar to 
#' \link[epimutacions]{epimutations} 
#' with the particularity that when 
#' is more than one case sample,  
#' the remaining case samples are included as controls. 
#' @param methy a GenomicRatioSet object
#'  containing the samples for the 
#' analysis. See the constructor function 
#' \link[minfi]{GenomicRatioSet}, 
#' \link[minfi]{makeGenomicRatioSetFromMatrix}. 
#' @param method a character string naming the 
#' outlier detection method to be used. 
#' This can be set as: 
#' \code{"manova"}, \code{"mlm"}, 
#' \code{"isoforest"}, \code{"mahdistmcd"}, 
#' \code{"barbosa"} and \code{beta}. 
#' The default is \code{"manova"}. 
#' For more information see \strong{Details}. 
#' @param epi_params the parameters for each method. 
#' See the function \link[epimutacions]{epi_parameters}.
#' @param BPPARAM (\code{"BiocParallelParam"}) 
#' \link[BiocParallel]{BiocParallelParam} object
#' to configure parallelization execution. 
#' By default, execution is non-parallel. 
#' @param verbose logical. If TRUE additional 
#' details about the procedure will provide to the user. 
#' The default is TRUE.
#' @param ... Further parameters passed to `epimutations`
#' @details The function compares a case 
#' sample against a control panel 
#' to identify epimutations in the given 
#' sample. First, the DMRs are identified using the 
#' \link[bumphunter]{bumphunter} approach. 
#' After that, CpGs in those DMRs are 
#' tested in order to detect regions
#' with CpGs being outliers.  
#' For that,  different anomaly 
#' detection methods can be selected:  
#'  * Multivariate Analysis of Variance 
#'  (\code{"manova"}). \link[stats]{manova}
#'  * Multivariate Linear Model (\code{"mlm"})
#'  * Isolation Forest (\code{"isoforest"}) 
#'  \link[isotree]{isolation.forest}
#'  * Robust Mahalanobis Distance (\code{"mahdistmcd"}) 
#'  \link[robustbase]{covMcd}
#'  * Barbosa (\code{"barbosa"})
#' @return The function returns an object of class tibble 
#' containing the outliers regions.  
#' The results are composed by the following columns: 
#' * \code{epi_id}: the name of the anomaly detection method that 
#' has been used to detect the epimutation
#' * \code{sample}: the name of the sample where the epimutation was found.
#' * \code{chromosome}, \code{start} and \code{end}: 
#' indicate the location of the epimutation.
#' * \code{sz}: the number of base pairs in the region.
#' * \code{cpg_n}: number of CpGs in the region.
#' * \code{cpg_ids}: differentially methylated CpGs names.
#' * \code{outlier_score}: 
#'    * For method \code{manova} it provides the approximation 
#'    to F-test and the Pillai score, separated by \code{/}.
#'    * For method \code{mlm} it provides the approximation to 
#'    F-test and the R2 of the model, separated by \code{/}.
#'    * For method \code{isoforest} it provides the 
#'    magnitude of the outlier score.
#'    * For methods \code{barbosa} and \code{mahdistmcd} is filled with NA.
#' * \code{outlier_significance}: 
#'    * For methods \code{manova}, \code{mlm}, and \code{isoforest} 
#'    it provides the p-value obtained from the model.
#'    * For method \code{barbosa} and \code{mahdistmcd} is filled with NA.
#' * \code{outlier_direction}: indicates the direction of 
#' the outlier with \code{"hypomethylation"} and \code{"hypermethylation"}
#'    * For \code{manova}, \code{mlm}, \code{isoforest}, and \code{mahdistmcd} 
#'    it is computed from the values obtained from bumphunter.
#'    * For \code{barbosa} it is computed from the location of 
#'    the sample in the reference distribution (left vs. right outlier).
#' @examples 
#' data(GRset)
#' \donttest{
#' epimutations_one_leave_out(GRset[,c(1:5,11)], method = "manova")
#' }
#' @importFrom methods is
#' @importFrom BiocParallel SerialParam bplapply 
#' @export 
epimutations_one_leave_out <- function(methy, method = "manova", 
                                       epi_params = epi_parameters(), 
                                       BPPARAM = BiocParallel::SerialParam(),
                                       verbose = TRUE, ...){
  
  if(!is(methy, "GenomicRatioSet"))
  {
    stop("Input data 'methy' must be a 'GenomicRatioSet'. 
         'makeGenomicRatioSetFromMatrix' function from 'minfi' package 
         can be useful to create a 'GenomicRatioSet' class object")
  }
  
  rst <- do.call(rbind, 
                 BiocParallel::bplapply(colnames(methy), 
                                        function(case){
    case_samples <- methy[, case]
    control_panel <-  methy[, !colnames(methy) %in% case]
    epimutations(case_samples, 
                 control_panel, 
                 method, 
                 epi_params = epi_params,
                 verbose = verbose, ...)
  }, BPPARAM = BPPARAM))
  return(rst)
}

