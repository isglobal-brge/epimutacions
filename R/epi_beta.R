#' @title Identifies epimutations based on a beta distribution.
#' 
#' @description \code{epi_beta} method models the DNA methylation data 
#' using a beta distribution. 
#' First, the beta distribution parameters of the 
#' reference population are precomputed and passed to the method. 
#' Then, we compute the probability of observing the methylation values 
#' of the case from the reference beta distribution. 
#' CpGs with p-values smaller than a threshold 
#' \code{pvalue_threshold} and with 
#' a methylation difference with the 
#' mean reference methylation
#' higher than \code{diff_threshold} 
#' are defined as outlier CpGs. Finally, 
#' epimutations are defined as a 
#' group of contiguous outlier CpGs.
#' @param beta_params matrix with the parameters of the reference 
#' beta distributions for each CpG in the dataset.
#' @param betas_case matrix with the methylation values for a case.
#' @param beta_mean beta values mean. 
#' @param betas a matrix containing the beta values for all samples.
#' @param case case sample name.  
#' @param controls control samples names.
#' @param annot annotation of the CpGs.
#' @param pvalue_threshold minimum p-value to consider a CpG an outlier.
#' @param diff_threshold minimum methylation difference 
#' between the CpG and the mean methylation to
#' consider a position an outlier. 
#' @param min_cpgs minimum number of CpGs to consider an epimutation.
#' @param maxGap maximum distance between two contiguous CpGs to 
#' combine them into an epimutation.
#' @return The function returns a data frame with 
#' the candidate regions to be epimutations.
#' 
#' @importFrom stats pbeta
#' 
epi_beta <-  function(beta_params, beta_mean, betas_case, case, controls, 
                    betas, annot, pvalue_threshold, diff_threshold, 
                    min_cpgs = 3, maxGap)
{

    ## Compute p-value for case
    if (requireNamespace("purrr", quietly = TRUE)) {
        pvals <- purrr::pmap_dbl(list(betas_case, beta_params[, 1],
                                        beta_params[, 2]),
                                function(x, shape1, shape2)
                                    stats::pbeta(x, shape1, shape2))
    } else {
        stop("'purrr' package not avaibale")
    }
    
    names(pvals) <- rownames(betas_case)
    
    ## Select CpGs with difference in mean methylation higher than threshold
    diff_vec <- abs(betas_case - beta_mean) > diff_threshold
    pvals <- pvals[diff_vec]
    
    ## Remove NAs
    pvals <- pvals[!is.na(pvals)]
    
    ## Hypomethylation regions
    negCpGs <- pvals[pvals < pvalue_threshold]
    negGR <- annot[names(negCpGs)]
    negGR$pvals <- negCpGs
    negRegs <- defineRegions(negGR, case, controls, betas, maxGap, up = FALSE)
    
    
    ## Hypermethylation regions
    posCpGs <- 1 - pvals[1 - pvals < pvalue_threshold]
    posGR <- annot[names(posCpGs)]
    posGR$pvals <- posCpGs
    posRegs <- defineRegions(posGR, case, controls, betas, maxGap)
    
    df <- rbind(posRegs, negRegs)
    
    if (nrow(df) > 0) {
        df <- subset(df, cpg_n >= min_cpgs, drop = FALSE)
    }
    df
}



#' @title  Model methylation as a beta distribution
#' @param x Matrix of methylation expressed as a beta. 
#' CpGs are in columns and samples in rows.   
#' @return Beta distribution. 
#' 
#' @importFrom matrixStats colVars
#' @importFrom bumphunter clusterMaker
#' @importFrom BiocGenerics start
#' @importFrom BiocGenerics start end width
#'                                       
getBetaParams <- function(x)
{
    xbar <- colMeans(x, na.rm = TRUE)
    s2 <- matrixStats::colVars(x, na.rm = TRUE)
    term <- (xbar * (1 - xbar)) / s2
    alpha.hat <- xbar * (term - 1)
    beta.hat <- (1 - xbar) * (term - 1)
    return(cbind(alpha.hat, beta.hat))
} 



defineRegions <- function(regGR, case, controls, betas, maxGap, up = TRUE)
{

    if (requireNamespace("GenomeInfoDb", quietly = TRUE)) {
        regGR <- sort(regGR)
        cl <- bumphunter::clusterMaker(GenomeInfoDb::seqnames(regGR),
                                    BiocGenerics::start(regGR),
                                    maxGap = maxGap)
        reg_list <- lapply(unique(cl), function(i) {
            cpgGR <- regGR[cl == i]
            rang <- range(cpgGR)
            data.frame( chromosome = as.character(GenomeInfoDb::seqnames(rang)),
                        start = BiocGenerics::start(rang),
                        end = BiocGenerics::end(rang),
                        sz = BiocGenerics::width(rang),
                        cpg_n = length(cpgGR),
                        cpg_ids = paste(names(cpgGR), collapse = ",", sep = ""),
                        outlier_score = mean(cpgGR$pvals),
                        outlier_direction = ifelse(up, "hypermethylation",
                                                        "hypomethylation"),
                        pvalue = NA,
                        adj_pvalue =  NA,
                        delta_beta = abs(mean(betas[names(cpgGR), controls]) -
                                mean(betas[names(cpgGR), case])),
                        sample = NA )
            })
    } else {
        stop("'GenomeInfoDb' package not avaibale")
    }
    
    if (length(reg_list) == 0) {
        df <- data.frame( chromosome = character(), start = numeric(),
                            end = numeric(), sz = numeric(),
                            cpg_n = numeric(), cpg_ids = character(),
                            outlier_score = numeric(),
                            outlier_direction = character(),
                            pvalue = numeric(), adj_pvalue = numeric(),
                            delta_beta = numeric() )
    } else {
        Reduce(rbind, reg_list)
    }
} 