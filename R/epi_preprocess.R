#' @export
#' @title Preprocess methylation array
#' @description The \code{epi_preprocess} function reads 
#' Illumina methylation sample 
#' sheet for case samples and it merges them with 
#' \link[minfi]{RGChannelSet} reference panel.
#'  The final dataset is normalized using minfi package preprocess methods. 
#' @param cases_dir the base directory from which the search is started. 
#' @param reference_panel an \link[minfi]{RGChannelSet} 
#' object containing the reference 
#' panel (controls) samples.
#' @param pattern What pattern is used to identify a sample sheet file. 
#' @param normalize a character string 
#' specifying the selected preprocess method. 
#' For more information see \strong{Details} or 
#' [minfi package user's Guide](shorturl.at/nwIN2). 
#' It can be set as: \code{"raw"}, \code{"illumina"}, 
#' \code{"swan"}, \code{"quantile"},
#' \code{"noob"} or \code{"funnorm"}.)
#' @param norm_param the parameters for each preprocessing method. 
#' See the function \link[epimutacions]{norm_parameters}.
#' @param verbose logical. If TRUE additional 
#' details about the procedure will provide to the user. 
#' The default is FALSE. 
#' @details 
#'  The \code{epi_preprocess} function reads Illumina methylation sample 
#' sheet for case samples and it merges them with 
#' \link[minfi]{RGChannelSet} reference panel.
#'  The final dataset is normalized using 
#'  different minfi package preprocess methods: 
#' * \code{"raw"}: \link[minfi]{preprocessRaw}
#' * \code{"illumina"}: \link[minfi]{preprocessIllumina}
#' * \code{"swan"}: \link[minfi]{preprocessSWAN}
#' * \code{"quantile"}: \link[minfi]{preprocessQuantile}
#' * \code{"noob"}: \link[minfi]{preprocessNoob}
#' * \code{"funnorm"}: \link[minfi]{preprocessFunnorm}
#' @return \code{epi_preprocess} function returns a 
#' \link[minfi]{GenomicRatioSet} object
#' containing case and control (reference panel) samples.  
#' @examples 
#'
#' # The reference panel for this example is available in 
#' #epimutacionsData (ExperimentHub) package
#' library(ExperimentHub)
#' eh <- ExperimentHub()
#' query(eh, c("epimutacionsData"))
#' reference_panel <- eh[["EH6691"]]
#' cases_dir <- system.file("extdata", package = "epimutacionsData")
#' #Preprocessing
#' \dontrun{
#' epi_preprocess(cases_dir, 
#'                reference_panel, 
#'                pattern = "SampleSheet.csv")
#'                }
#' 
#' 
#' 
#' @importFrom minfi read.metharray.sheet read.metharray.exp combineArrays 
#' preprocessRaw preprocessIllumina preprocessSWAN preprocessQuantile
#' preprocessNoob preprocessFunnorm mapToGenome ratioConvert

epi_preprocess <-function(cases_dir, 
                          reference_panel, 
                          pattern = "csv$",  
                          normalize = "raw", 
                          norm_param = norm_parameters(), 
                          verbose = FALSE){
  
  if(is.null(cases_dir))
  {
    stop("The argument 'cases_dir' must be introduced")
    
  }
  
  if (!requireNamespace("methods")) 
    stop("'methods' package not available")
  
  avail <- c("raw", "illumina", 
             "swan", "quantile", 
             "noob", "funnorm")
  normalize <- charmatch(normalize, avail)
  normalize <- avail[normalize]
  if(is.na(normalize)) 
    stop("Invalid normalisation ('normalize') method was selected'")
  
  #Reading case samples idat files
  targets <- minfi::read.metharray.sheet(base = cases_dir, 
                                         pattern = pattern)
  if(is.null(targets)){
    warning("There is not any sample sheet in the base directory")
    RGset_cases <- minfi::read.metharray.exp(base = cases_dir)
  }else{
    RGset_cases <- minfi::read.metharray.exp(targets = targets)
  }
  
 #Reference panel

  if(!is(reference_panel, "RGChannelSet")){
    stop("Reference panel must be a 'RGChannelSet' class object")
  }
 

 #Merge reference panel and case samples
 RGset <- minfi::combineArrays(reference_panel, RGset_cases,
                 outType = c("IlluminaHumanMethylation450k",
                           "IlluminaHumanMethylationEPIC"),
                 verbose = TRUE)
 
 #Preprocesing 
 c("raw","illumina", "swan", "quantile", "noob", "funnorm")
 if(normalize == "raw"){
   Mset <- minfi::preprocessRaw(RGset)
 }else if(normalize == "illumina"){
   Mset <-minfi::preprocessIllumina(RGset, 
                                    bg.correct = 
                                      norm_param$illumina$bg.correct,
                                    normalize = 
                                      norm_param$illumina$normalize,
                                    reference = 
                                      norm_param$illumina$reference)
 }else if(normalize == "swan"){
   Mset <- minfi::preprocessSWAN(RGset, 
                                 verbose = verbose)
 }else if(normalize == "quantile"){
   GRset <- minfi::preprocessQuantile(
              RGset, 
              fixOutliers = 
                norm_param$quantile$fixOutliers,
              removeBadSamples = 
                norm_param$quantile$removeBadSamples,
              badSampleCutoff = 
                norm_param$quantile$badSampleCutoff,
              quantileNormalize = 
                norm_param$quantile$quantileNormalize,
              stratified = 
                norm_param$quantile$stratified,
              mergeManifest = 
                norm_param$quantile$mergeManifest, 
              sex = 
                norm_param$quantile$sex,
              verbose = verbose)
   
 }else if(normalize == "noob"){
   Mset <- minfi::preprocessNoob(RGset, 
                                 offset = norm_param$noob$offset,
                                 dyeCorr = norm_param$noob$dyeCorr, 
                                 verbose = verbose,
                                 dyeMethod = norm_param$noob$dyeMethod)
 }else if(normalize == "funnorm"){
   GRset <- minfi::preprocessFunnorm(RGset, 
                                     nPCs = norm_param$funnorm$nPCs,
                                     sex = norm_param$funnorm$sex, 
                                     bgCorr = norm_param$funnorm$bgCorr,
                                     dyeCorr = norm_param$funnorm$dyeCorr,
                                     keepCN = norm_param$funnorm$keepCN,
                                     ratioConvert = TRUE,
                                     verbose = verbose)
 }
 
 #Create GenomicRatioSet object
 if(normalize %in% c("raw",
                     "illumina", 
                     "swan", 
                     "noob")){
   GMset <- minfi::mapToGenome(Mset, 
                               mergeManifest = FALSE)
   GRset <- minfi::ratioConvert(GMset, 
                                what = "beta", 
                                keepCN = FALSE)
   
 }
 return(GRset)
}

