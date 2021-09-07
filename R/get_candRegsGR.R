get_candRegsGR <- function(){
  eh <- ExperimentHub()
  query(eh, c("epimutacionsData"))
  return(eh[["EH6692"]])
}