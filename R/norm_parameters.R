norm_parameters <- function(illumina = list("bg.correct" = TRUE, "normalize" = c("controls", "no"), reference = 1),
                            quantile = list("fixOutliers" = TRUE, "removeBadSamples" = FALSE,
                                            "badSampleCutoff" = 10.5, "quantileNormalize" = TRUE,
                                            "stratified" = TRUE, "mergeManifest" = FALSE, "sex" = NULL), 
                            noob = list("offset" = 15, "dyeCorr" = TRUE, "dyeMethod" = c("single", "reference")), 
                            funnorm = list("nPCs" = 2, "sex" = NULL, "bgCorr" = TRUE,
                                           "dyeCorr" = TRUE, "keepCN" = FALSE)){
                              
                              
                              return(list("illumina" = illumina ,"quantile" = quantile, "noob" = noob, "funnorm" = funnorm))
                            }