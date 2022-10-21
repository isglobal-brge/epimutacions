##' Non-parametric, Asymptotic P-values for Multivariate Linear Models
##' 
##' Fits a multivariate linear model and computes 
##' test statistics and asymptotic 
##' P-values for predictors in a non-parametric manner. 
##' 
##' A \code{Y} matrix is obtained after transforming 
##' (optionally) and centering 
##' the original response variables.
##'  Then, the multivariate fit obtained by 
##' \code{\link{lm}} can be used to 
##' compute sums of squares, pseudo-F statistics and asymptotic 
##' P-values for the terms specified by the \code{formula} 
##' in a non-parametric manner. 
##' @param formula object of class "\code{\link{formula}}" 
##' (or one that can be 
##' coerced to that class): a symbolic description 
##' of the model to be fitted. 
##' @param data an optional data frame, 
##' list or environment (or object coercible 
##' by \code{\link{as.data.frame}} to a data frame) 
##' containing the variables in 
##' the model. If not found in data, 
##' the variables are taken from 
##' \code{environment(formula)}, 
##' typically the environment from which \code{mlm} 
##' is called.
##' @param transform transformation 
##' of the response variables: "\code{none}", 
##' "\code{sqrt}" or "\code{log}". 
##' Default is "\code{none}".
##' @param contrasts an optional list.
##'  See \code{contrasts.arg} in 
##' \code{\link{model.matrix.default}}. 
##' Default is "\code{\link{contr.sum}}" 
##' for ordered factors and "\code{\link{contr.poly}}" 
##' for unordered factors. 
##' Note that this is different from the default setting in 
##' \code{\link{options}("contrasts")}.
##' @param subset subset of predictors for which 
##' summary statistics will be 
##' reported. Note that this is different 
##' from the "\code{subset}" 
##' argument in \code{\link{lm}}.
##' @param fit logical. If \code{TRUE} 
##' the multivariate fit on transformed and 
##' centered responses is returned.
##' 
##' @return \code{mlm} returns an object 
##' of \code{\link{class}} \code{"MLM"}, 
##' a list containing:
##' \item{call}{the matched call.}
##' \item{aov.tab}{ANOVA table with Df, Sum Sq, Mean Sq, F values, 
##' partial R-squared and P-values.}
##' \item{precision}{the precision in P-value computation.}
##' \item{transform}{the transformation 
##' applied to the response variables.}
##' \item{na.omit}{incomplete cases removed 
##' (see \code{\link{na.omit}}).}
##' \item{fit}{if \code{fit = TRUE} the multivariate 
##' fit done on the transformed 
##' and centered response variables is also returned.}
##' 
##' @seealso \code{\link{lm}}, \code{\link[car]{Anova}}
##' 
##' @author Diego Garrido-Martín
##' 
##' @importFrom stats model.frame model.response lm.fit
##' .getXlevels
##' 
##' 
##' 
mlm <- function(formula, data, transform = "none", 
                contrasts = NULL, subset = NULL, fit = FALSE)
{
    ## Checks
    transform <- match.arg(transform, c("none", "sqrt", "log"))
    
    ## Save call, build model frame, obtain responses
    cl <- match.call()
    m <- match(c("formula", "data"), names(cl), 0L)
    mf <- cl[c(1L, m)]
    mf$drop.unused.levels <- TRUE
    mf$na.action <- "na.omit"
    # The rows with at least one NA either in Y or X
    # (only considering variables used in the formula)
    # will be removed before transforming/centering
    mf[[1L]] <- quote(stats::model.frame)
    mf <- eval(mf, parent.frame())
    mt <- attr(mf, "terms")
    response <- stats::model.response(mf, "numeric")
    
    ## Checks
    if (NCOL(response) < 2) {
        stop("The number of response variables should be >= 2")
    }
    if (length(attr(mt, "term.labels")) < 1) {
        stop("The model should contain at least one predictor
                        (excluding the intercept)")
    }
    
    ## Transform and center responses, update model frame
    if (transform == "none") {
        Y <- response
    } else if (transform == "sqrt") {
        if (any(response < 0)) {
            stop("'sqrt' transformation requires all response values >= 0")
        }
        Y <- sqrt(response)
    } else if (transform == "log") {
        if (any(response <= 0)) {
            stop("'log' transformation requires all response values > 0")
        }
        Y <- log(response)
    }
    Y <- scale(Y, center = TRUE, scale = FALSE)
    mf[[1L]] <- Y
    
    ## Define contrasts
    if (is.null(contrasts)) {
        contrasts <- list(unordered = "contr.sum", ordered = "contr.poly")
        dc <- attr(mt, "dataClasses")[-1]
        contr.list <- lapply(
            dc,
            FUN = function(k) {
                # No contrast for quantitative predictors
                # Sum contrasts for unordered categorical predictors
                # Polynomial contrasts for ordered categorical predictors
                contr.type <- switch(k, "factor" = contrasts$unordered,
                                    "ordered" = contrasts$ordered)
                return(contr.type)
            }
        )
        contr.list <- contr.list[!unlist(lapply(contr.list, is.null))]
    } else {
        contr.list <- contrasts
    }
    
    ## Build model matrix
    X <- model.matrix(mt, mf, contr.list)
    
    ## Fit lm
    lmfit <- stats::lm.fit(X, Y)
    class(lmfit) <- c("mlm", "lm")
    lmfit$na.action <- attr(mf, "na.action")
    lmfit$contrasts <- attr(X, "contrasts")
    lmfit$xlevels <- stats::.getXlevels(mt, mf)
    lmfit$call <- cl[c(1L, m)]
    lmfit$call[[1L]] <- quote(lm)
    if (length(contr.list) > 0)
        lmfit$call$contrasts <- quote(contr.list)
    lmfit$terms <- mt
    lmfit$model <- mf
    
    ## Compute sums of squares, df's, pseudo-F
    ## statistics, partial R2s and eigenvalues
    stats <- mlmtst(fit = lmfit, X = X, subset = subset)
    SS <- stats$SS
    df <- stats$df
    f.tilde <- stats$f.tilde
    r2 <- stats$r2
    e <- stats$e
    
    ## Compute P-values
    l <- length(df) # SS[l], df[l] correspond to Residuals
    pv.acc <- mapply(p.asympt, ss = SS[-l], df = df[-l], 
                    MoreArgs = list(lambda = e))
    
    ## ANOVA table
    stats.l <- list(df, SS, SS / df, f.tilde, r2, pv.acc[1,])
    cmat <- do.call(rbind, stats.l)
    cmat <- t(cmat)
    colnames(cmat) <- c("Df", "Sum Sq", "Mean Sq", "F value", "R2", "Pr(>F)")
    
    ## Output
    out <- list( "call" = cl, "aov.tab" = cmat, "type" = type,
                "precision" = pv.acc[2,], "transform" = transform,
                "na.omit" = lmfit$na.action )
    if (fit) {
        out$fit <- lmfit
    }
    
    ## Update class
    class(out) <- c('MLM', class(out))
    return(out)
}

