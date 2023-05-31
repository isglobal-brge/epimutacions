integrateReferenceQuantile <- function(betas, beta_params){
  
  joint <- cbind(beta_params, betas)
  corrected <- apply(joint, 1, function(x) normalizeBeta(x[-c(1:2)], x[1:2]))
  corrected <- t(corrected)
  dimnames(corrected) <- dimnames(betas)
  corrected
}


normalizeBeta <- function(vals, beta_params){
  
  ### Select values from quantiles 5-95%
  center <- vals > quantile(vals, 0.05, na.rm = TRUE) & vals < quantile(vals, 0.95, na.rm = TRUE)
  valsf <- vals[center]
  bpars <- getBetaParams(matrix(vals, ncol = 1))
  
  ps <- pbeta(vals, bpars[1], bpars[2])
  vals_pit <- qbeta(ps, beta_params[1], beta_params[2])
  
}