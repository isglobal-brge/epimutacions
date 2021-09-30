#'  Non-parametric, Asymptotic P-values for Multivariate Linear Models
#' 
#'  Fits a multivariate linear model and computes test statistics and asymptotic 
#'  P-values for predictors in a non-parametric manner. 
#' 
#'  A \code{Y} matrix is obtained after transforming (optionally) and centering 
#'  the original response variables. Then, the multivariate fit obtained by 
#'  \code{\link{lm}} can be used to compute sums of squares (type-I, type-II or 
#'  type-III), pseudo-F statistics and asymptotic P-values for the terms specified
#'  by the \code{formula} in a non-parametric manner. The designations "type-II" 
#'  and "type-III" correspond exactly to those used in \code{\link[car]{Anova}}. 
#'  "type-I" refers to sequential sums of squares.
#' 
#'  @param formula object of class "\code{\link{formula}}" (or one that can be 
#'  coerced to that class): a symbolic description of the model to be fitted. 
#'  @param data an optional data frame, list or environment (or object coercible 
#'  by \code{\link{as.data.frame}} to a data frame) containing the variables in 
#'  the model. If not found in data, the variables are taken from 
#'  \code{environment(formula)}, typically the environment from which \code{mlm} 
#'  is called.
#'  @param transform transformation of the response variables: "\code{none}", 
#'  "\code{sqrt}" or "\code{log}". Default is "\code{none}".
#'  @param type type of sum of squares: "\code{I}", "\code{II}" or "\code{III}". 
#'  Default is "\code{II}".
#'  @param contrasts an optional list. See \code{contrasts.arg} in 
#'  \code{\link{model.matrix.default}}. Default is "\code{\link{contr.sum}}" 
#'  for ordered factors and "\code{\link{contr.poly}}" for unordered factors. 
#'  Note that this is different from the default setting in \code{\link{options}("contrasts")}.
#'  @param subset subset of predictors for which summary statistics will be 
#'  reported. Note that this is different from the "\code{subset}" argument in \code{\link{lm}}.
#'  @param fit logical. If \code{TRUE} the multivariate fit on transformed and 
#'  centered responses is returned.
#' 
#'  @details  ORIGINAL AUTHOR: Diego Garrido-Martín
#'  The original mlm from mlm v0.8.3
#'  at GitHub: https://github.com/dgarrimar/mlm
#' 
#'  @return \code{mlm} returns an object of \code{\link{class}} "MLM", a list containing:
#'  \item{call}{the matched call.}
#'  \item{aov.tab}{ANOVA table with Df, Sum Sq, Mean Sq, F values, 
#'  partial R-squared and P-values.}
#'  \item{type}{the type of sum of squares (\code{"I"}, \code{"II"} or \code{"III"}).}
#'  \item{precision}{the precision in P-value computation.}
#'  \item{transform}{the transformation applied to the response variables.}
#'  \item{na.omit}{incomplete cases removed (see \code{\link{na.omit}}).}
#'  \item{fit}{if \code{fit = TRUE} the multivariate fit done on the transformed 
#'  and centered response variables is also returned.}
#' 
#' @seealso \code{\link{lm}}, \code{\link[car]{Anova}}
#' 
#' @author Diego Garrido-Martín
#' 
#' @import stats
#' 
#' @export
#'  @details 
#'
#' 
mlm <- function(formula, data, transform = "none", type = "II", 
                contrasts = NULL, subset = NULL, fit = FALSE){
  
  ## Checks
  transform <- match.arg(transform, c("none", "sqrt", "log"))
  type <- match.arg(type, c("I", "II", "III"))
  
  ## Save call, build model frame, obtain responses
  cl <- match.call()
  m <- match(c("formula", "data"), names(cl), 0L)
  mf <- cl[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf$na.action = "na.omit"
  # The rows with at least one NA either in Y or X 
  # (only considering variables used in the formula) 
  # will be removed before transforming/centering
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())               
  mt <- attr(mf, "terms")
  response <- model.response(mf, "numeric")
  
  ## Checks
  if (NCOL(response) < 2){
    stop("The number of response variables should be >= 2")
  }
  if (length(attr(mt, "term.labels")) < 1) {
    stop("The model should contain at least one predictor (excluding the intercept)")
  }
  
  ## Transform and center responses, update model frame
  if(transform == "none"){
    Y <- response
  } else if (transform == "sqrt"){
    if (any(response < 0)) {
      stop("'sqrt' transformation requires all response values >= 0")
    }
    Y <- sqrt(response)
  } else if (transform == "log"){
    if (any(response <= 0)) {
      stop("'log' transformation requires all response values > 0")
    }
    Y <- log(response)
  }
  Y <- scale(Y, center = TRUE, scale = FALSE)
  mf[[1L]] <- Y
  
  ## Define contrasts
  if(is.null(contrasts)){
    contrasts <- list(unordered = "contr.sum", ordered = "contr.poly")
    dc <- attr(mt, "dataClasses")[-1]
    contr.list <- lapply(dc, FUN = function(k){
      # No contrast for quantitative predictors
      # Sum contrasts for unordered categorical predictors
      # Polynomial contrasts for ordered categorical predictors
      contr.type <- switch(k, "factor" = contrasts$unordered,
                           "ordered" = contrasts$ordered)
      return(contr.type)
    })
    contr.list <- contr.list[!unlist(lapply(contr.list, is.null))]
  } else {
    contr.list <- contrasts
  }
  
  ## Build model matrix
  X <- model.matrix(mt, mf, contr.list)
  
  ## Fit lm
  lmfit <- lm.fit(X, Y)
  class(lmfit) <- c("mlm", "lm")
  lmfit$na.action <- attr(mf, "na.action")
  lmfit$contrasts <- attr(X, "contrasts")
  lmfit$xlevels <- .getXlevels(mt, mf)
  lmfit$call <- cl[c(1L, m)]
  lmfit$call[[1L]] <- quote(lm)
  if(length(contr.list) > 0) lmfit$call$contrasts <- quote(contr.list)
  lmfit$terms <- mt
  lmfit$model <- mf
  
  ## Compute sums of squares, df's, pseudo-F statistics, partial R2s and eigenvalues 
  stats <- mlmtst(fit = lmfit, X = X, type = type, subset = subset)
  SS <- stats$SS
  df <- stats$df
  f.tilde <- stats$f.tilde
  r2 <- stats$r2
  e <- stats$e
  
  ## Compute P-values
  l <- length(df) # SS[l], df[l] correspond to Residuals 
  pv.acc <- mapply(p.asympt, ss = SS[-l], df = df[-l], MoreArgs = list(lambda = e))  
  
  ## ANOVA table
  stats.l <- list(df, SS, SS/df, f.tilde, r2, pv.acc[1, ])
  cmat <- data.frame()
  for(i in seq(along = stats.l)) {
    for(j in names(stats.l[[i]])){
      cmat[j, i] <- stats.l[[i]][j]
    }
  }
  cmat <- as.matrix(cmat)
  colnames(cmat) <- c("Df", "Sum Sq", "Mean Sq", "F value", "R2", "Pr(>F)")
  
  ## Output
  out <- list("call" = cl,
              "aov.tab" = cmat,
              "type" = type,
              "precision" = pv.acc[2, ],
              "transform" = transform,
              "na.omit" = lmfit$na.action)
  if(fit){
    out$fit <- lmfit
  }
  
  ## Update class
  class(out) <- c('MLM', class(out))
  return(out)
}

#' Asymptotic P-values
#' 
#' Computes asymptotic P-values given the numerator of the pseudo-F statistic, 
#' its degrees of freedom and the eigenvalues of the residual covariance matrix.
#' 
#' @param ss numerator of the pseudo-F statistic.
#' @param lambda eigenvalues of the residual covariance matrix.
#' @param df degrees of freedom of the numerator of the pseudo-F statistic.
#' @param eps the desired level of accuracy.
#' @param eps.updt factor by which \code{eps} is updated to retry execution of
#' algorithm AS 204 when it fails with fault indicator 4, 5 or 9.
#' @param eps.stop if \code{eps > eps.stop}, execution of algorithm AS 204 is 
#' not retried and the function raises an error. Default is \code{1e-10}.
#' 
#' @details  ORIGINAL AUTHOR: Diego Garrido-Martín
#' The original mlm from mlm v0.8.3
#' at GitHub: https://github.com/dgarrimar/mlm
#' 
#' @return A vector containing the P-value and the level of accuracy. 
#' 
#' @seealso \code{\link{AS204}}
#' 
#' @author Diego Garrido-Martín
#' 
#' @keywords internal
#' 
#' 
p.asympt <- function(ss, df, lambda, eps = 1e-14, eps.updt = 2, eps.stop = 1e-10){
  
  pv  <- AS204(c = ss, lambda = lambda, mult = rep(df, length(lambda)), eps = eps)
  while (is.null(pv)) {
    eps <- eps * eps.updt
    if(eps > eps.stop){
      stop(sprintf("Precision of asymptotic P-value > %.2e", eps.stop))
    }
    pv  <- AS204(c = ss, lambda = lambda, mult = rep(df, length(lambda)), eps = eps)
  }
  if(pv < eps){
    pv <- eps
  }
  return(c(pv, eps))
}

#' Sums of Squares and Pseudo-F Statistics from a Multivariate Fit
#' 
#' Computes the sum of squares, degrees of freedom, pseudo-F statistics and
#' partial R-squared for each predictor from a multivariate \code{fit}. 
#' It also returns the eigenvalues of the residual covariance matrix.
#' 
#' Different types of sums of squares (i.e. "\code{I}", "\code{II}" and 
#' "\code{III}") are available.
#' 
#' @param fit multivariate fit obtained by \code{\link{lm}}.
#' @param X design matrix obtained by \code{\link{model.matrix}}. 
#' @param type type of sum of squares ("\code{I}", "\code{II}" or "\code{III}"). 
#' Default is "\code{II}".
#' @param subset subset of predictors for which summary statistics will be reported.
#' Note that this is different from the "\code{subset}" argument in \code{\link{lm}}.
#' @param tol \code{e[e/sum(e) > tol]}, where \code{e} is the vector of eigenvalues
#' of the residual covariance matrix. Required to prevent long running times of 
#' algorithm AS 204. Default is 0.001 to ensure minimal loss of accuracy.
#' 
#' @details  ORIGINAL AUTHOR: Diego Garrido-Martín
#' The original mlm from mlm v0.8.3
#' at GitHub: https://github.com/dgarrimar/mlm
#' 
#' @return A list containing:
#' \item{SS}{sums of squares for all predictors (and residuals).}
#' \item{df}{degrees of freedom for all predictors (and residuals).}
#' \item{f.tilde}{pseudo-F statistics for all predictors.}
#' \item{r2}{partial R-squared for all predictors.}
#' \item{e}{eigenvalues of the residual covariance matrix.}
#'
#' @seealso \code{\link{AS204}}
#' 
#' @author Diego Garrido-Martín
#'
#' @keywords internal
#' 
#' 
mlmtst <- function(fit, X, type = "II", subset = NULL, tol = 1e-3){
  
  ## Residual sum-of-squares and cross-products (SSCP) matrix
  SSCP.e <- crossprod(fit$residuals) 
  
  ## Residual sum-of-squares and df
  SS.e <- sum(diag(SSCP.e))
  df.e <- fit$df.residual # df.e <- (n-1) - sum(df)
  
  ## Total sum-of-squares
  SS.t <- sum(diag(crossprod(fit$model[[1L]])))
  
  ## Partial sums-of-squares
  terms <- attr(fit$terms, "term.labels") # Model terms
  n.terms <- length(terms)
  
  if(!is.null(subset)) {
    if(all(subset %in% terms)) {
      iterms <- which(terms %in% subset)
    }else {
      stop(sprintf("Unknown terms in subset: %s",
                   paste0("'", subset[which(! subset %in% terms)], "'", 
                          collapse = ", ")))
    }
  } else {
    iterms <- 1:n.terms
  }
  asgn <- fit$assign
  
  df <- SS <- numeric(n.terms) # Initialize empty
  names(df) <- names(SS) <- terms
  
  if (type == "I"){
    
    effects <- as.matrix(fit$effects)[seq_along(asgn), , drop = FALSE]
    
    for (i in iterms) {
      subs <- which(asgn == i) 
      SS[i] <- sum(diag(crossprod(effects[subs, , drop = FALSE])))
      df[i] <- length(subs)
    }
    
  } else {
    
    sscp <- function(L, B, V){
      LB <- L %*% B
      crossprod(LB, solve(L %*% tcrossprod(V, L), LB))
    }
    
    B <- fit$coefficients     # Coefficients
    V <- solve(crossprod(X))  # V = (X'X)^{-1}
    p <- nrow(B)
    I.p <- diag(p)
    
    # In contrast to car::Anova, intercept 
    # information is not returned for
    # type III sums-of-squares
    
    if (type == "III"){
      
      for (i in iterms){
        subs <- which(asgn == i) 
        L <- I.p[subs, , drop = FALSE] # Hypothesis matrix
        SS[i] <- sum(diag(sscp(L, B, V)))
        df[i] <- length(subs)
      }
      
    } else {
      
      is.relative <- function(term1, term2, factors) {
        all( !( factors[, term1] & ( !factors[, term2] ) ) )
      }
      
      fac <- attr(fit$terms, "factors") 
      for (i in iterms){
        term <- terms[i]
        subs.term <- which(asgn == i)
        if(n.terms > 1) { # Obtain relatives
          relatives <- (1:n.terms)[-i][sapply(terms[-i], 
                                              function(term2) 
                                                is.relative(term, term2, fac))]
        } else { 
          relatives <- NULL
        }
        subs.relatives <- NULL
        for (relative in relatives){
          subs.relatives <- c(subs.relatives, which(asgn == relative))
        }
        L1 <- I.p[subs.relatives, , drop = FALSE] # Hyp. matrix (relatives) 
        if (length(subs.relatives) == 0) {
          SSCP1 <- 0
        } else {
          SSCP1 <- sscp(L1, B, V)
        }
        L2 <- I.p[c(subs.relatives, subs.term), , drop = FALSE] # Hyp. matrix (relatives + term) 
        SSCP2 <- sscp(L2, B, V)
        SS[i] <- sum(diag(SSCP2 - SSCP1))
        df[i] <- length(subs.term)
      }
    }
  }
  
  ## subset
  if(!is.null(subset)){
    SS <- SS[iterms]
    df <- df[iterms]
  }
  
  ## pseudo-F
  f.tilde <- SS/SS.e*df.e/df
  
  ## r.squared
  R2 <- (SS.t - SS.e)/SS.t
  # R2adj <- 1-( (1-R2)*(n-1) / df.e ) 
  r2 <- SS/SS.t
  # r2adj <- 1-( (1-r2)*(n-1) / df.e )
  
  # Get eigenvalues from cov(R)*(n-1)/df.e
  e <- eigen(SSCP.e/df.e, symmetric = T, only.values = T)$values
  e <- e[e/sum(e) > tol]
  
  return(list("SS" = c(SS, "Residuals" = SS.e),
              "df" = c(df, "Residuals" = df.e),
              "f.tilde" = f.tilde, "r2" = r2, "e" = e))
}

#' Algorithm AS 204
#' 
#' Distribution of a positive linear combination of \eqn{\chi^2} random variables.
#' 
#' Algorithm AS 204 evaluates the expression
#' \deqn{
#' P [X < c] = P [ \sum_{j=1}^n \lambda_j \chi^2(m_j, \delta^2_j) < c ]
#' }
#' where \eqn{\lambda_j} and \eqn{c} are positive constants and 
#' \eqn{\chi^2(m_j, \delta^2_j)} represents an independent \eqn{\chi^2} 
#' random variable with \eqn{m_j} degrees of freedom and non-centrality
#' parameter \eqn{\delta^2_j}. This can be approximated by the truncated series 
#' \deqn{
#' \sum_{k=0}^{K-1} a_k P [\chi^2(m+2k) < c/\beta]
#' } 
#' where \eqn{m = \sum_{j=1}^n m_j} and \eqn{\beta} is an arbitrary constant 
#' (as given by argument "mode"). 
#' 
#' The \code{C++} implementation of algorithm AS 204 used here is identical 
#' to the one employed by the \code{\link[CompQuadForm]{farebrother}} method
#' in the \code{CompQuadForm} package, with minor modifications.
#' 
#' @param c value point at which distribution is to be evaluated.
#' @param lambda the weights \eqn{\lambda_j}.
#' @param mult the multiplicities \eqn{m_j}.
#' @param delta the non-centrality parameters \eqn{\delta^2_j}.
#' @param maxit the maximum number of terms \eqn{K} (see Details).
#' @param eps the desired level of accuracy.
#' @param mode if "\code{mode}" > 0 then \eqn{\beta=mode\lambda_{min}},
#' otherwise \eqn{\beta=2/(1/\lambda_{min}+1/\lambda_{max})}.
#' 
#' @details  ORIGINAL AUTHOR: Diego Garrido-Martín
#' The original mlm from mlm v0.8.3
#' at GitHub: https://github.com/dgarrimar/mlm
#' 
#' @return The function returns the probability \eqn{P[X > c] = 1 - P[X < c]} 
#' if the AS 204 fault indicator is 0 (see Note below), and \code{NULL} if 
#' the fault indicator is 4, 5 or 9, as the corresponding faults can be
#' corrected by increasing "\code{eps}". Other faults raise an error. 
#' 
#' @note The algorithm AS 204 defines the following fault indicators:
#' \strong{-j)} one or more of the constraints \eqn{\lambda_j > 0}, 
#' \eqn{m_j > 0} and \eqn{\delta^2_j \ge 0} is not satisfied.
#' \strong{1)} non-fatal underflow of \eqn{a_0}.
#' \strong{2)} one or more of the constraints \eqn{n > 0}, 
#' \eqn{c > 0}, \eqn{maxit > 0} and \eqn{eps > 0} is not satisfied.
#' \strong{3)} the current estimate of the probability is < -1.
#' \strong{4)} the required accuracy could not be obtained in \eqn{maxit} iterations.
#' \strong{5)} the value returned by the procedure does not satisfy 
#' \eqn{0 \le P [X < c] \le 1}.
#' \strong{6)} the density of the linear form is negative.
#' \strong{9)} faults 4 and 5.
#' \strong{10)} faults 4 and 6. 
#' \strong{0)} otherwise.
#' 
#' 
#' @references 
#' P. Duchesne, P. Lafaye de Micheaux, Computing the distribution of quadratic forms: 
#' Further comparisons between the Liu-Tang-Zhang approximation and exact methods, 
#' Computational Statistics and Data Analysis, Vol. 54, (2010), 858-862 
#' 
#' Farebrother R.W., Algorithm AS 204: The distribution of a Positive Linear Combination 
#' of chi-squared random variables, Journal of the Royal Statistical Society, 
#' Series C (applied Statistics), Vol. 33, No. 3 (1984), 332-339
#' 
#' @seealso \link[CompQuadForm]{farebrother}
#' 
#' @author Diego Garrido-Martín
#' 
#' @useDynLib mlm ruben
#' 
#' @keywords internal
#' 
AS204 <- function (c, lambda, mult = rep(1, length(lambda)), delta = rep(0, length(lambda)),
                   maxit = 100000, eps = 1e-14, mode = 1) {
  
  out <- .C("ruben", lambda = as.double(lambda), mult = as.integer(mult), 
            delta = as.double(delta), n = as.integer(length(lambda)), 
            c = as.double(c), mode = as.double(mode), maxit = as.integer(maxit), 
            eps = as.double(eps), dnsty = as.double(0), ifault = as.integer(0), 
            res = as.double(0), PACKAGE = "mlm")
  
  if(out$ifault == 0){
    return(1 - out$res)
  } else if (out$ifault %in% c(4, 5, 9)){
    return(NULL)
  } else {
    stop(sprintf("Algorithm AS 204 failed with fault indicator: %s", out$ifault))
  }
}

