#' @export
#' @title Settings for parameters of \code{epimutations} 
#' and \code{epimutations_one_leave_out} functions
#' @description  Allow the user to set the values 
#' of the parameters to compute the functions 
#' \link[epimutacions]{epimutations} 
#' and \link[epimutacions]{epimutations_one_leave_out}.
#' @param manova,mlm,iForest,mahdist,quantile,beta method selected 
#' in the function \link[epimutacions]{epimutations}. 
#' @param pvalue_cutoff the threshold p value to select 
#' which CpG regions are outliers in \code{manova}, 
#' \code{mlm} and \code{beta} methods.
#' @param outlier_score_cutoff The outlier score 
#' threshold to identify outliers CpGs in
#' isolation forest (\code{iForest}) method. Default is \code{0.5}. 
#' @param ntrees number of binary trees to build for the model build by 
#' isolation forest (\code{iForest}) method. 
#' Default is \code{100}. 
#' @param nsamp the number of subsets used for initial 
#' estimates in the Minimum Covariance Determinant 
#' which is used to compute the Robust Mahalanobis 
#' distance (\code{mahdist}). 
#' It can be set as:
#' \code{"best"}, \code{"exact"}, or \code{"deterministic"}. 
#' For \code{nsamp = "best"} exhaustive enumeration is done, 
#' as long as the number of trials does not exceed 100'000. 
#' For \code{nsamp = "exact"} exhaustive enumeration 
#' will be attempted however many samples are needed. 
#' In this case, a warning message may be displayed 
#' saying that the computation can take a very long time.
#' For \code{nsamp = "deterministic"}. 
#' For more information see \link[robustbase]{covMcd}.
#' Default is \code{"deterministic"}. 
#' @param window_sz the maximum distance between 
#' CpGs to be considered in the same DMR. 
#' This parameter is used in \code{quantile} (default: 1000). 
#' @param qsup,qinf,offset_abs The upper and lower quantiles (threshold) 
#' to consider a CpG an outlier when using \code{quantile} method, 
#' as well as the offset to consider (defaults: 0.005, 0.995, 0.15).
#' @param diff_threshold Minimum methylation difference between 
#' the CpG and the mean methylation to consider a position an outlier. 
#' @details Invoking \code{epi_parameters()} with 
#' no arguments returns return a list with the default values. 
#' @return  the function returns a list of all 
#' set parameters for each method used in 
#' \link[epimutacions]{epimutations} and 
#' \link[epimutacions]{epimutations_one_leave_out} functions.
#' @examples 
#'    #Default set of parameters
#'    epi_parameters()
#'    #change p value for manova method
#'    epi_parameters(manova = list("pvalue_cutoff" = 0.01))
#' @export

epi_parameters <- function( manova = list("pvalue_cutoff" = 0.05), 
                            mlm = list("pvalue_cutoff" = 0.05),
                            iForest = list("outlier_score_cutoff" = 0.7,
                                            "ntrees" = 100),
                            mahdist = list("nsamp" = "deterministic"),
                            quantile = list( "window_sz"=1000,
                                            "offset_abs"=0.15,
                                            "qsup"=0.995,
                                            "qinf"=0.005 ),
                            beta =  list("pvalue_cutoff" = 1e-6,
                                        "diff_threshold" = 0.1)) 
{
    return( list( "manova" = manova ,
                    "mlm" = mlm,
                    "iForest" = iForest,
                    "mahdist" = mahdist,
                    "quantile" = quantile,
                    "beta" = beta 
                )
            )
}