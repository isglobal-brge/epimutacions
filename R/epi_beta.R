epi_beta <-  function(beta_params, betas_case, annot, pvalue_threshold, 
                      min_cpgs = 3, maxGap){
 
  
  ## Compute p-value for case
  mat <- cbind(betas_case, t(beta_params))
  pvals <- apply(mat, 1, function(x) {
    p <- pbeta(x[1], x[2], x[3], lower.tail = TRUE)
    p <- min(p, 1 - p)
  })

  ## Select Significant CpGs
  sigCpGs <- names(pvals[!is.na(pvals) & pvals < pvalue_threshold])
  sigGR <- annot[sigCpGs]
  overs <- findOverlaps(sigGR, candRegsGR)
  regs <- table(to(overs))
  selRegs <- names(regs[regs >= min_cpgs])
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


## Model methylation as a beta distribution
getBetaParams <- function(x){
  xbar <- colMeans(x)
  s2 <- colVars(x)
  term <- (xbar*(1-xbar))/s2
  alpha.hat <- xbar*(term-1)
  beta.hat <- (1-xbar)*(term-1)
  return(cbind(alpha.hat, beta.hat))
} 

getBetaParams <- function(mat){
  
  ## Add offset to avoid undefined values in optimization
  mat[mat == 0] <- 0.001
  mat[mat == 1] <- 0.999

  ## Get beta distribution parameters
  beta_params <- getBetaParams(mat)
  beta_params
}



# reticulate::source_python("/home/carlos/Repositories/combined-pvalues/cpv/pipeline.py")
# pipeline("pvalue", "None", "None", 333, "out", "None", 0.05, "None", "~/Repositories/combined-pvalues/examples/file.bed")
