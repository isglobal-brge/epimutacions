#' @export
#' @title Settings for parameters of epi_preprocess function 
#' @description  \code{norm_parameters} function allows the user 
#' to set the values of the parameters to compute the functions
#' \code{epi_preprocess}.
#' @param illumina,quantile,noob,funnorm preprocess method 
#' selected in the function \link[epimutacions]{epi_preprocess}. 
#' @param bg.correct logical. If TRUE background 
#' correction will be performed in \code{"illumina"} method. Default is TRUE.
#' @param normalize logical. If TRUE control 
#' normalization will be performed in 
#' \code{"illumina"} method. 
#' @param reference numeric. 
#' The reference array for control normalization in 
#' \code{"illumina"} method. 
#' @param fixOutliers logical. If TRUE low outlier 
#' Meth and Unmeth signals will be fixed
#' in \code{"quantile"} method. Default is TRUE. 
#' @param removeBadSamples logical. 
#' If TRUE bad samples will be removed.  
#' @param badSampleCutoff a numeric specifying 
#' the cutoff to label samples
#' as 'bad' in \code{"quantile"} method. Default is 10.5.
#' @param quantileNormalize logical. If TRUE quantile 
#' normalization will be performed in 
#' \code{"quantile"} method. Default is TRUE.
#' @param stratified logical. 
#' If TRUE quantile normalization will be performed
#' within region strata in \code{"quantile"} method. 
#' Default is TRUE. 
#' @param mergeManifest logical. If TRUE the 
#' information in the associated manifest
#' package will be merged into the output 
#' object in \code{"quantile"} method. 
#' Default is FALSE. 
#' @param sex an optional numeric vector containing 
#' the sex of the samples in \code{"quantile"} method. 
#' @param offset a numeric specifying an offset for 
#' the normexp background correction
#' in \code{"noob"} method. Default is 15. 
#' @param dyeCorr logial. Dye correction will 
#' be done in \code{"noob"} and \code{"funnorm"} methods. 
#' Default is TRUE. 
#' @param dyeMethod specify the dye 
#' bias correction to be done, single sample 
#' approach or a reference array in \code{"noob"} method. 
#' @param nPCs numeric specifying 
#' the number of principal components 
#' from the control probes PCA in 
#' \code{"funnorm"} method. Default is 2. 
#' @param sex an optional numeric vector 
#' containing the sex of the samples in 
#' \code{"quantile"} and \code{"funnorm"} methods.
#' @param bgCorr logical. 
#' If TRUE NOOB background correction will be done prior 
#' to functional normalization. 
#' in \code{"funnorm"} method. Default is TRUE. 
#' @param keepCN logical. If TRUE copy number estimates 
#' will be kept in \code{"funnorm"} method. 
#' Default is FALSE.
#' @details Invoking \code{epi_parameters()} with no 
#' arguments returns a list with the
#' default values for each normalization parameter. 
#' @return  the function returns a list of all 
#' set parameters for each normalization method used in 
#' \code{epi_peprocess}. 
#' @examples 
#'  #Default set of parameters
#'  norm_parameters()
#'  #change p value for manova method
#'  norm_parameters(illumina = list("bg.correct" = FALSE))
#'  
norm_parameters <- function(illumina = 
                              list("bg.correct" = TRUE, 
                                   "normalize" = c("controls", "no"), 
                                   reference = 1),
                            quantile = list("fixOutliers" = TRUE, 
                                            "removeBadSamples" = FALSE,
                                            "badSampleCutoff" = 10.5, 
                                            "quantileNormalize" = TRUE,
                                            "stratified" = TRUE, 
                                            "mergeManifest" = FALSE, 
                                            "sex" = NULL), 
                            noob = list("offset" = 15, 
                                        "dyeCorr" = TRUE, 
                                        "dyeMethod" = c("single", 
                                                        "reference")), 
                            funnorm = list("nPCs" = 2, 
                                           "sex" = NULL, 
                                           "bgCorr" = TRUE,
                                           "dyeCorr" = TRUE, 
                                           "keepCN" = FALSE)){
                              
                              
                              return(list("illumina" = illumina,
                                          "quantile" = quantile,
                                          "noob" = noob,
                                          "funnorm" = funnorm))
                            }