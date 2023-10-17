#' #' @title  PCA correction
#' #' 
#' #' @description  Compute PCA correction
#' #' 
#' #' @param cases a GenomicRatioSet object containing the case samples.
#' #' See the constructor function \link[minfi]{GenomicRatioSet}, 
#' #' \link[minfi]{makeGenomicRatioSetFromMatrix}. 
#' #' @param controls a GenomicRatioSet object containing the 
#' #' control panel (control panel).
#' #' @param prange probe.range
#' #' 
#' #' @return A list with cases and controls corrected by PCs
#' #' 
#' #' 
#' #' @examples 
#' #' 
#' #' data(res.epi.manova)
#' #' #Annotate the epimutations
#' #' 
#' #' 
#' PCA_correction <- function(cases, controls, prange = 40000) {
#'     
#'     ## Compute residuals of pcs
#'     p <- intersect( featureNames(cases), featureNames(controls))
#'     ## found 51 overlapping featureNames then combine the sets
#'     all_samples <- Biobase::combine( cases[p, ], controls[p, ])
#'     
#'     beta <- meffil:::impute.matrix(minfi::getBeta(all_samples), margin = 1)
#'     ndim <- isva::EstDimRMT(beta, FALSE)$dim + 1
#'     pcs <- meffil::meffil.methylation.pcs(minfi::getBeta(all_samples), probe.range = prange)
#'     m <- getM(all_samples)
#'     res <- residuals(limma::lmFit(m, pcs[, seq_len(ndim)]), m)
#'     beta <- ilogit2(res)
#'     assay(all_samples) <- beta
#'     
#'     gr_cases <- all_samples[ , colnames(cases)]
#'     gr_controls <- all_samples[ , colnames(controls)]
#'     
#'     return(list(cases = gr_cases, 
#'                 controls = gr_controls))
#'     
#' }