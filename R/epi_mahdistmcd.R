#' @title Identifies epimutations using Robust Mahalanobis distance
#' @description  This function identifies regions with CpGs being outliers
#' using the Minimum Covariance Determinant (MCD) estimator 
#' (\link[robustbase]{covMcd}) to compute the Mahalanobis distance. 
#' @param mixture beta values matrix. Samples in columns and
#' CpGs in rows.
#' @param nsamp the number of subsets used for initial estimates in the MCD. 
#' It can be set as:
#' \code{"best"}, \code{"exact"}, or \code{"deterministic"}. 
#' @details The implementation of the method here is based 
#' on the discussion in this
#' thread of [Cross Validated](https://cutt.ly/Kka0M87)
#' @return The function returns the computed Robust Mahalanobis distance.
#' 

epi_mahdistmcd <- function(mixture, 
                           nsamp = c("best", "exact", "deterministic")) {
	nsamp <- charmatch(nsamp, c("best", "exact", "deterministic"))
	nsamp <- c("best", "exact", "deterministic")[nsamp]
	if(is.na(nsamp)) {
		stop("Argument 'nsamp' shuld be 'best', 'exact', 'deterministic'")
	}
	
	# Transpose input methylation beta matrix to have the samples as rows and 
	# the CpGs as columns
	mixture <- t(mixture)
	
	# Run the MCD model on the transposed matrix
	mcd <- robustbase::covMcd(mixture, nsamp = nsamp)
	
	# Get MCD stimate of location
	mean_mcd <- mcd$center
	
	# Get MCD estimate scatter
	cov_mcd <- mcd$cov
	
	# Get inverse of scatter
	cov_mcd_inv <- solve(cov_mcd)
	
	# Compute the robust distance between samples
	robust_dist <- apply(mixture, 1, function(x){
		x <- (x - mean_mcd)
		dist <- sqrt((t(x)  %*% cov_mcd_inv %*% x))
		return(dist)
	})
	
	return(data.frame(ID = names(robust_dist), statistic = robust_dist))
}

#' @title  Creates a data frame containing the results 
#' obtained from Robust Mahalanobis distance
#' @description Creates a data frame containing the
#' genomic regions, statistics and direction for the DMRs.
#' @param bump a DMR obtained from \link[bumphunter]{bumphunter}
#' (i.e. a row from \link[bumphunter]{bumphunter} method result).
#' @param outliers  the robust distance computed by 
#'  \link[epimutacions]{epi_mahdistmcd} function results. 
#' @param case a character string specifying the case sample name. 
#' @returns The function returns a data frame containing 
#' the following information for each DMR: 
#' * genomic ranges
#' * DMR base pairs
#' * number and name of CpGs in DMR
#' * statistics: 
#'     * Outlier score
#'     * Outlier significance
#'     * Outlier direction
#'  * Sample name
#' 
#' For more information about the output see \link[epimutacions]{epimutations}.

res_mahdistmcd <- function(case, bump, outliers) {
	bump$outlier <- case %in% outliers
	bump$outlier_score <- NA
	bump$pvalue <- NA 
	bump$adj_pvalue <- NA
	bump$outlier_direction <- NA
	bump <- bump[bump$outlier, ]
	bump <- bump[ , c("chromosome", 
	                  "start", 
	                  "end", "sz", "cpg_n", "cpg_ids", "outlier_score",
	                  "outlier_direction", "pvalue", "adj_pvalue",  "sample")]	
	return(bump)
}

