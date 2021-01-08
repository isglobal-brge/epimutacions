#' @export
epi_params <- function(manova = list("pvalue_cutoff" = 0.05), mlm = list("pvalue_cutoff" = 0.05), 
                       isoforest = list("outlier_score_cutoff" = 0.5), mahdistmcd = list("nsamp" = "deterministic"),
                       barbosa = list("window_sz" = 10, "offset_mean" = 0.15, "offset_abs" = 0.1), qn = list("window_sz" = 10, "qn_th" = 3)){
                      
  return(list("manova" = manova ,"mlm" = mlm,"isoforest" = isoforest, "mahdistmcd" = mahdistmcd,"barbosa" = barbosa,"qn" = qn))
}






  

