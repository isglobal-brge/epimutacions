integrateReferenceQuantile <- function(betas, beta_params){
  
  joint <- cbind(beta_params, betas)
  corrected <- apply(joint, 1, function(x) normalizeBeta(x[-c(1:2)], x[1:2]))
  corrected <- t(corrected)
  dimnames(corrected) <- dimnames(betas)
  corrected
}


normalizeBeta <- function(vals, beta_params){
  
  ### Select values from quantiles 5-95%
  iqr <- stats::IQR(vals, na.rm = TRUE)
  outs <- vals < quantile(vals, 0.25, na.rm = TRUE) - 1.5*iqr | vals > quantile(vals, 0.75, na.rm = TRUE) + 1.5*iqr
  outs[is.na(outs)] <- TRUE
  valsf <- vals[!outs]
  bpars <- getBetaParams(matrix(valsf, ncol = 1))
  
  mean_ori <- beta_params[1] / (beta_params[1] + beta_params[2])
  mean_new <- bpars[1] / (bpars[1] + bpars[2])
  
  ps <- pbeta(vals, bpars[1], bpars[2])
  vals_pit <- qbeta(ps, beta_params[1], beta_params[2])
  
  vals_pit[outs] <- vals[outs] + (mean_ori - mean_new)
  vals_pit
}