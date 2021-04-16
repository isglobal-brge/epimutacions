epi_preprocess <-function(cases_dir, reference_panel,  normalise = "raw", norm_param = norm_parameters(), verbose = FALSE){
  
  if(is.null(cases_dir))
  {
    stop("The argument 'cases_dir' must be introduced")
    
  }
  
  avail <- c("raw","illumina", "swan", "quantile", "noob", "funnorm")
  normalise <- charmatch(normalise, avail)
  normalise <- avail[normalise]
  if(is.na(normalise)) stop("Invalid normalisation ('normalise') method was selected'")
  
  #Reading case samples idat files
  targets <- minfi::read.metharray.sheet(cases_dir)
  if(is.null(targets)){
    warning("There is not any sample sheet in the base directory")
    RGset_cases <- minfi::read.metharray.exp(cases_dir)
  }
  RGset_cases <- minfi::read.metharray.exp(targets = targets)
  
 #Reference panel

  if(class(reference_panel) != "RGChannelSet"){
    stop("Reference panel must be a 'RGChannelSet' class object")
  }
 

 #Merge reference panel and case samples
 RGset <- minfi::combineArrays(reference_panel, RGset_cases,
                 outType = c("IlluminaHumanMethylation450k",
                           "IlluminaHumanMethylationEPIC"),
                 verbose = TRUE)
 
 #Preprocesing 
 c("raw","illumina", "swan", "quantile", "noob", "funnorm")
 if(normalise == "raw"){
   Mset <- minfi::preprocessRaw(RGset)
 }else if(normalise == "illumina"){
   Mset <- minfi::preprocessIllumina(RGset, bg.correct = norm_param$illumina$bg.correct,
                                     normalize = norm_param$illumina$normalize,
                                     reference = norm_param$illumina$reference)
 }else if(normalise == "swan"){
   Mset <- minfi::preprocessSWAN(RGset, verbose = verbose)
 }else if(normalise == "quantile"){
   GRset <- minfi::preprocessQuantile(RGset, fixOutliers = norm_param$quantile$fixOutliers,
                                      removeBadSamples = norm_param$quantile$removeBadSamples,
                                      badSampleCutoff = norm_param$quantile$badSampleCutoff,
                                      quantileNormalize = norm_param$quantile$quantileNormalize,
                                      stratified = norm_param$quantile$stratified,
                                      mergeManifest = norm_param$quantile$mergeManifest, 
                                      sex = norm_param$quantile$sex,
                                      verbose = verbose)
   
 }else if(normalise == "noob"){
   Mset <- minfi::preprocessNoob(RGset, offset = norm_param$noob$offset,
                                 dyeCorr = norm_param$noob$dyeCorr, verbose = verbose,
                                 dyeMethod = norm_param$noob$dyeMethod)
 }else if(normalise == "funnorm"){
   GRset <- minfi::preprocessFunnorm(RGset, nPCs = norm_param$funnorm$nPCs,
                                     sex = norm_param$funnorm$sex, 
                                     bgCorr = norm_param$funnorm$bgCorr,
                                     dyeCorr = norm_param$funnorm$dyeCorr,
                                     keepCN = norm_param$funnorm$keepCN,
                                     ratioConvert = TRUE,
                                     verbose = verbose)
 }
 
 #Create GenomicRatioSet object
 if(normalise %in% c("raw","illumina", "swan", "noob")){
   GMset <- minfi::mapToGenome(Mset, mergeManifest = FALSE)
   GRset <- minfi::ratioConvert(GMset, what = "beta", keepCN = FALSE)
   
 }
 return(GRset)
}

