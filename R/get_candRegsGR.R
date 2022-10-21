#' @title Candidate regions to be epimutations 
#' @description  Load candidate regions to be epimutations 
#' from \code{epimutacionsData} package in \code{ExperimentHub}. 
#' @return  The function returns a GRanges object containing 
#' the candidate regions.
#' 
#' @import epimutacionsData
#' @importFrom ExperimentHub ExperimentHub

get_candRegsGR <- function(){
    if (!requireNamespace("ExperimentHub"))
        stop("'ExperimentHub' package not available")
    if (!requireNamespace("AnnotationHub"))
        stop("'AnnotationHub' package not available")
    eh <- ExperimentHub()
    query(eh, c("epimutacionsData"))
    return(eh[["EH6692"]])
}