#' @title Identifies epimutations based on a beta distribution.
#' 
#' @description This method models the DNA methylation data using a beta distribution. First, 
#' the beta distribution parameters of the reference population are precomputed and passed
#' to the method. Then, we compute the probability of observing the methylation values 
#' of the case from the reference beta distribution. CpGs with p-values smaller than a threshold 
#' \code{pvalue_threshold} and with a methylation difference with the mean reference methylation
#' higher than \code{diff_threshold} are defined as outlier CpGs. Finally, epimutations are defined as group
#' of contiguous outlier CpGs.
#' 
#' @param beta_params Matrix with the parameters of the reference beta distributions for 
#' each CpG in the dataset.
#' @param betas_case Matrix with the methylation values for a case
#' @param annot Annotation of the CpGs
#' @param pvalue_threshold Minimum p-value to consider a CpG an outlier
#' @param diff_threshold Minimum methylation difference between the CpG and the mean methylation to
#' consider a position an outlier. 
#' @param min_cpgs Minimum number of CpGs to consider an epimutation
#' @param maxGap Maximum distance between two contiguous CpGs to combine them into an epimutation
#' @return The function returns a data frame with the regions candidates to be
#' epimutations.
epi_beta <-  function(beta_params, beta_mean, betas_case, annot, pvalue_threshold, 
                      diff_threshold, min_cpgs = 3, maxGap){
  
  
  ## Compute p-value for case
  pvals <- purrr::pmap_dbl(list(betas_case, beta_params[, 1], beta_params[, 2]), 
                           function(x, shape1, shape2) pbeta(x, shape1, shape2))
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
  negRegs <- defineRegions(negGR, maxGap, up = FALSE)
  
  
  ## Hypermethylation regions
  posCpGs <- 1 - pvals[1 - pvals < pvalue_threshold]
  posGR <- annot[names(posCpGs)]
  posGR$pvals <- posCpGs
  posRegs <- defineRegions(posGR, maxGap)
  
  df <- rbind(posRegs, negRegs)
  
  if (nrow(df) > 0){
    df <- subset(df, N_CpGs >= min_cpgs, drop = FALSE )
  }
  df
}



#' @title  Model methylation as a beta distribution
#' @param x Matrix of methylation expressed as a beta. CpGs are in columns and samples in rows.                       
getBetaParams <- function(x){
  xbar <- colMeans(x, na.rm = TRUE)
  s2 <- matrixStats::colVars(x, na.rm = TRUE)
  term <- (xbar*(1-xbar))/s2
  alpha.hat <- xbar*(term-1)
  beta.hat <- (1-xbar)*(term-1)
  return(cbind(alpha.hat, beta.hat))
} 

defineRegions <- function(regGR, maxGap, up = TRUE){
  
  regGR <- sort(regGR)
  cl <- bumphunter::clusterMaker(seqnames(regGR), start(regGR), maxGap = maxGap)
  reg_list <- lapply(unique(cl), function(i){
    cpgGR <- regGR[cl == i]
    rang <- range(cpgGR)
    data.frame(chromosome = seqnames(rang), start = start(rang), 
               end = end(rang),
               length = width(rang), N_CpGs = length(cpgGR), 
               cpg_ids = paste(names(cpgGR), collapse = ",", sep = ""),
               outlier_score = mean(cpgGR$pvals),
               outlier_significance = NA,
               outlier_direction = ifelse(up, "hypermethylation", "hypomethylation"),
               sample = NA
    )
  })
  if (length(reg_list) == 0){
    df <-   data.frame(chromosome = character(), start = numeric(), 
                       end = numeric(),
                       length = numeric(), N_CpGs = numeric(), 
                       cpg_ids = character(),
                       outlier_score = numeric(),
                       outlier_significance = numeric(),
                       outlier_direction = character()
    )
  } else {
    Reduce(rbind, reg_list)
  }
} 