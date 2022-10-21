##' Sums of Squares and Pseudo-F Statistics from a Multivariate Fit
##' 
##' Computes the sum of squares, 
##' degrees of freedom, pseudo-F statistics and
##' partial R-squared for 
##' each predictor from a multivariate \code{fit}. 
##' It also returns the eigenvalues of 
##' the residual covariance matrix.
##' 
##' Different types of sums of squares 
##' (i.e. "\code{I}", "\code{II}" and 
##' "\code{III}") are available.
##' 
##' @param fit multivariate fit obtained by \code{\link{lm}}.
##' @param X design matrix obtained by \code{\link{model.matrix}}. 
##' @param subset subset of predictors for which summary 
##' statistics will be reported. 
##' Note that this is different from the 
##' "\code{subset}" argument in \code{\link{lm}}.
##' @param tol \code{e[e/sum(e) > tol]}, where \code{e} 
##' is the vector of eigenvalues
##' of the residual covariance matrix. 
##' Required to prevent long running times of 
##' algorithm AS 204. Default is 0.001 to 
##' ensure minimal loss of accuracy.
##' 
##' @return A list containing:
##' \item{SS}{sums of squares for all predictors (and residuals).}
##' \item{df}{degrees of freedom for all predictors (and residuals).}
##' \item{f.tilde}{pseudo-F statistics for all predictors.}
##' \item{r2}{partial R-squared for all predictors.}
##' \item{e}{eigenvalues of the residual covariance matrix.}
##'
##' @seealso \code{\link{AS204}}
##' 
##' @author Diego Garrido-Martín
##'
##' @keywords internal
##' 
mlmtst <- function(fit, X,  subset = NULL,  tol = 1e-3)
{
    
    ## Residual sum-of-squares and cross-products (SSCP) matrix
    SSCP.e <- crossprod(fit$residuals)
    
    ## Residual sum-of-squares and df
    SS.e <- sum(diag(SSCP.e))
    df.e <- fit$df.residual
    
    ## Total sum-of-squares
    SS.t <- sum(diag(crossprod(fit$model[[1L]])))
    
    ## Partial sums-of-squares
    terms <- attr(fit$terms, "term.labels") # Model terms
    n.terms <- length(terms)
    
    if (!is.null(subset)) {
        if (all(subset %in% terms)) {
            iterms <- which(terms %in% subset)
        } else {
            stop(sprintf( "Unknown terms in subset: %s",
                paste0("'", subset[which(!subset %in% terms)], "'", 
                        collapse = ", ")))
        }
    } else {
        iterms <- seq_len(n.terms)
    }
    asgn <- fit$assign
    
    sscp <- function(L, B, V) {
        LB <- L %*% B
        crossprod(LB, solve(L %*% tcrossprod(V, L), LB))
    }
    
    B <- fit$coefficients     # Coefficients
    V <- solve(crossprod(X))
    p <- nrow(B)
    I.p <- diag(p)
    
    term <- terms[1]
    subs.term <- which(asgn == 1)
    relatives <- NULL
    subs.relatives <- NULL
    L1 <- I.p[subs.relatives, , drop = FALSE]
    SSCP1 <- 0
    L2 <- I.p[c(subs.relatives, subs.term), , drop = FALSE]
    # Hyp. matrix (relatives + term)
    SSCP2 <- sscp(L2, B, V)
    SS <- sum(diag(SSCP2 - SSCP1))
    df <- length(subs.term)
    
    ## pseudo-F
    f.tilde <- SS / SS.e * df.e / df
    
    ## r.squared
    R2 <- (SS.t - SS.e) / SS.t
    r2 <- SS / SS.t
    
    e <- eigen(SSCP.e / df.e, symmetric = TRUE, only.values = TRUE)$values
    e <- e[e / sum(e) > tol]
    
    return(list( "SS" = c(SS, "Residuals" = SS.e),
                "df" = c(df, "Residuals" = df.e),
                "f.tilde" = f.tilde,
                "r2" = r2, 
                "e" = e ))
}