#' @title  Annotate the results of  epimutations or 
#' epimutations_one_leave_out functions
#' 
#' @description  Information about close genes and regulatory 
#' elements for epimutations.
#' 
#' @param epi_results a data frame object containing the output 
#' from \code{epimutations} or
#' \code{epimutations_one_leave_out} functions. 
#' @param db a character string containing  
#' the Illumina annotation package 
#' used to annotate the CpGs.
#' @param build a character string containing the genomic 
#' build where the epimutations are mapped. 
#' The default is GRCh37 (\code{build = "37"}). To use GRCh38 
#' set \code{built} to \code{NULL}. 
#' @param ... Further arguments passed to \code{annotate_cpg}.
#' 
#' @return The function returns the 
#' input object \code{epi_results}
#' with additional columns containing the  information about
#' the genes or overlapping regulatory features. 
#' 
#' See \link[epimutacions]{annotate_cpg} and 
#' \link[epimutacions]{add_ensemble_regulatory}
#' for an in-depth description of these variables.
#' 
#' @examples 
#' 
#' data(res.epi.manova)
#' #Annotate the epimutations
#' 
#' #anno_results <- annotate_epimutations(res.epi.manova)
#' 
#' 
#' @export
annotate_epimutations <- function(epi_results, 
                        db = "IlluminaHumanMethylationEPICanno.ilm10b2.hg19", 
                        build = "37", ...)
{
    
    if (!requireNamespace(db, quietly = TRUE)) {
        stop( db, " package not available")
    }
    
    if (!requireNamespace(db, quietly = TRUE)) {
        stop( "'Homo.sapiens' package not available")
    }
    ## Add gene mapping and CpG island information
    epi_results <- annotate_cpg(epi_results, db = db,  ...)
    
    ## Add ENSEMBL regulatory regions
    epi_results <- add_ensemble_regulatory(epi_results, build = build)
    
    epi_results
}
