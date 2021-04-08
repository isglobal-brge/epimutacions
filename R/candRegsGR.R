#' \code{GenomicRanges} class object
#'
#' @details Regions were epimutations can be defined. We cluster the CpGs that were closer
#' than 1000Kb to their closest CpG. As a result, we defined 40408 epimutations regions. Due to
#' the microarray design, all possible epimutations should be included in a distinct 
#' epimutation region. 
#' 
#' Next, we mapped epimutations regions to hg38. We selected those regions whose mapping
#' between hg19 and hg38 builds was 100%. Next, we overlapped epimutation regions with 
#' ENCODE cREs v13. 
#' 
#' Notice that epimutation regions are defined as the longest region that an epimutation
#' can have, considering that all the CpGs are outlier. Most epimutations regions are 
#' short, so the epimutations detected will be equivalent to the epimutation region. Nonetheless,
#' some epimutation regions are very long (> 5Kb). In these cases, the detected epimutation will
#' comprise only a short region of the epimutation region. Therefore, the overlap with 
#' cREs should be checked. 
#' 
#' @usage data("candRegsGR")
#' @return a \code{GenomicRanges} object.
#' @examples
#' data("candRegsGR")
#' candRegsGR
#'
"candRegsGR"
