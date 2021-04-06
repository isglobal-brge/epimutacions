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
  df <- Reduce(rbind, epi_list)
  subset(df, N_CpGs >= min_cpgs )
}


## Model methylation as a beta distribution
getBetaParamsVec <- function(vec){
  
  llhd2 <- function(x,p) {
    ans <- -sum(log(dbeta(x, p[1], p[2])))
    ans
  }
  
  vec <- vec[!is.na(vec)]
  pIni <- c(1, 1)
  
  param <- nlm(llhd2, x = vec, p = pIni)$estimate
  param
}

getBetaParams <- function(mat){
  
  ## Add offset to avoid undefined values in optimization
  mat[mat == 0] <- 0.001
  mat[mat == 1] <- 0.999

  ## Get beta distribution parameters
  beta_params <- apply(mat, 1, getBetaParamsVec)
  beta_params
}



# reticulate::source_python("/home/carlos/Repositories/combined-pvalues/cpv/pipeline.py")
# pipeline("pvalue", "None", "None", 333, "out", "None", 0.05, "None", "~/Repositories/combined-pvalues/examples/file.bed")
