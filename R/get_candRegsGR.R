#' @title Candidate regions to be epimutations 
#' @description  Load candidate regions to be epimutations 
#' from \code{epimutacionsData} package in \code{ExperimentHub}. 
#' @import ExperimentHubData
#' @import epimutacionsData
#' @return  The function returns a GRanges object containing 
#' the candidate regions. 
#' 

get_candRegsGR <- function(){
  eh <- ExperimentHub::ExperimentHub()
  AnnotationHub::query(eh, c("epimutacionsData"))
  return(eh[["EH6692"]])
}