#' @title Obtains bumps beta values 
#' @description The function obtains beta values 
#' corresponding to the CpGs into DMRs.
#' @param bump the result from \link[bumphunter]{bumphunter}.
#' @param fd a data frame containing the genomic ranges for each CpGs.
#' @param betas a matrix containing the beta 
#' values for all CpGs in each sample.
#' @return The function returns a data 
#' frame containing the beta values
#' for each sample and CpG into DMR.
betas_from_bump <- function(bump, fd, betas)
{
    cpgs <- rownames(fd[fd$seqnames == bump$chr & fd$start >= bump$start &
                        fd$start <= bump$end,])
    return(betas[cpgs,])
}
