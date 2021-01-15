#' @export
#' @title Settings for EpiMutations and epimutacions parameters
#' @Description Allow the user to set the values of the parameters to compute the functions
#' EpiMutations and epimutacions
#' @param pvalue_cutoff the threshold p value to select outliers DMR in \code{"manova"}
#' and \code{"mlm"} methods
#' @param  outlier_score_cutoff The outlier score threshold to identify outliers DMR in
#' isolation forest (\code{"isoforest"}) method
#' @param nsamp The number of subsets used for initial estimates. It can be set as:
#' \code{"best"}, \code{"exact"}, or \code{"deterministic"}. 
#'  For nsamp = "best" exhaustive enumeration is done, as long as the number of trials does not exceed 100'000 (= nLarge). For "exact", exhaustive enumeration will be attempted however many samples are needed. In this case a warning message may be displayed saying that the computation can take a very long time.
#'  For "deterministic", the deterministic MCD is computed; as proposed by Hubert et al. (2012) it starts from the h most central observations of six (deterministic) estimators.

epi_parameters <- function(manova = list("pvalue_cutoff" = 0.05), mlm = list("pvalue_cutoff" = 0.05), 
                       isoforest = list("outlier_score_cutoff" = 0.5), mahdistmcd = list("nsamp" = "deterministic"),
                       barbosa = list("window_sz" = 10, "offset_mean" = 0.15, "offset_abs" = 0.1), qn = list("window_sz" = 10, "qn_th" = 3)){
                      
  return(list("manova" = manova ,"mlm" = mlm,"isoforest" = isoforest, "mahdistmcd" = mahdistmcd,"barbosa" = barbosa,"qn" = qn))
}






  

