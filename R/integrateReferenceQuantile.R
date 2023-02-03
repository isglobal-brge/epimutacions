integrateReferenceQuantile <- function(betas, quantiles){
  
  joint <- cbind(quantiles[, -c(1:2)], betas)
  corrected <- apply(joint, 1, function(x) normalizeQuantile(x[-c(1:19)], x[1:19]))
  corrected <- t(corrected)
  dimnames(corrected) <- dimnames(betas)
  corrected
}


normalizeQuantile <- function(vals, quantiles){
  ### Select values from quantiles 5-95%
  center <- vals > quantile(vals, 0.05, na.rm = TRUE) & vals < quantile(vals, 0.95, na.rm = TRUE)
  valsf <- vals[center & !is.na(vals)]
  valsQ <- preprocessCore::normalize.quantiles.use.target(matrix(valsf, ncol = 1), quantiles)
  moddf <- list(Q = valsQ, f = valsf)
  
  diff <- valsf - valsQ
  mod <- .lm.fit(cbind(matrix(valsf), 1), valsQ)
  bpars <- getBetaParams(matrix(valsf, ncol = 1))
  ps <- pbeta(vals, bpars[1], bpars[2])
  quants <- pmax(ps, 1 - ps)
  vals_out <- (1 - quants)/0.5*(vals*mod$coefficients[1] + mod$coefficients[2]) + (quants - 0.5)/0.5*(vals - mean(diff))
  
  vals_out[vals_out < 1e-3] <- 1e-3
  vals_out[vals_out > 1-1e-3] <- 1-1e-3
  vals_out
  
}