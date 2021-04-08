#' @title Identifies epimutations based on a beta distribution.
#' 
#' @description This method models the DNA methylation data using a beta distribution. First, 
#' the beta distribution parameters of the reference population are precomputed and passed
#' to the method. Then, we compute the probability of observing the methylation values 
#' of the case from the reference beta distribution. CpGs with p-values smaller than a threshold 
#' \code{pvalue_threshold} are defined as outlier CpGs. Finally, epimutations are defined as group
#' of contiguous outlier CpGs.
#' 
#' @param beta_params Matrix with the parameters of the reference beta distributions for 
#' each CpG in the dataset.
#' @param betas_case Matrix with the methylation values for a case
#' @param annot Annotation of the CpGs
#' @param pvalue_threshold Minimum p-value to consider a CpG an outlier
#' @param min_cpgs Minimum number of CpGs to consider an epimutation
#' @param maxGap Maximum distance between two contiguous CpGs to combine them into an epimutation
#' @return The function returns a data frame with the regions candidates to be
#' epimutations.
epi_beta <-  function(beta_params, betas_case, annot, pvalue_threshold, 
                      min_cpgs = 3, maxGap){
 
  
  ## Compute p-value for case
  pvals <- purrr::pmap_dbl(list(betas_case, beta_params[, 1], beta_params[, 2]), 
                           function(x, shape1, shape2){ 
                             p <- pbeta(x, shape1, shape2)
                             p <- min(p, 1 - p)
                             })
    names(pvals) <- rownames(betas_case)
  
  ## Select Significant CpGs
  sigCpGs <- names(pvals[!is.na(pvals) & pvals < pvalue_threshold])
  sigGR <- annot[sigCpGs]
  
  ## Overlap significant CpGs with epimutations regions to speed up epimutations definition
  overs <- findOverlaps(sigGR, candRegsGR)
  regs <- table(to(overs))
  selRegs <- names(regs[regs >= min_cpgs])
  
  ## Refine epimutations definition
  epi_list <- lapply(selRegs, function(idx){
    ov <- overs[to(overs) == idx]
    regGR <- sort(sigGR[from(ov)])
    
    cl <- bumphunter::clusterMaker(seqnames(regGR), start(regGR), maxGap = maxGap)
    reg_list <- lapply(unique(cl), function(i){
      cpgGR <- regGR[cl == i]
      rang <- range(cpgGR)
      data.frame(chromosome = seqnames(rang), start = start(rang), 
                 end = end(rang),
                 length = width(rang), N_CpGs = length(cpgGR), 
                 cpg_ids = paste(names(cpgGR), collapse = ",", sep = ""),
                 outlier_score = mean(pvals[names(cpgGR)]),
                 outlier_significance = NA,
                 outlier_direction = NA
      )
    })
    Reduce(rbind, reg_list)
    

  })
  if (length(epi_list) == 0){
    df <-   data.frame(chromosome = character(), start = numeric(), 
                       end = numeric(),
                       length = numeric(), N_CpGs = numeric(), 
                       cpg_ids = character(),
                       outlier_score = numeric(),
                       outlier_significance = numeric(),
                       outlier_direction = character()
    )
  } else {
    df <- Reduce(rbind, epi_list)
  }
  df <- subset(df, N_CpGs >= min_cpgs )
  df
}


# Model methylation as a beta distribution
#' @param x Matrix of methylation expressed as a beta. CpGs are in columns and samples in rows.                       
getBetaParams <- function(x){
  xbar <- colMeans(x, na.rm = TRUE)
  s2 <- matrixStats::colVars(x, na.rm = TRUE)
  term <- (xbar*(1-xbar))/s2
  alpha.hat <- xbar*(term-1)
  beta.hat <- (1-xbar)*(term-1)
  return(cbind(alpha.hat, beta.hat))
} 

