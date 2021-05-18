#' @title Identifies epimutations using MANOVA
#' @description  This function identifies regions with CpGs being outliers
#' using \link[stats]{manova} approach. 
#' @param mixture beta values matrix. Samples in columns and
#' CpGs in rows.
#' @param  model design (or model) matrix.
#' @param case_id a character string specifying the name of the case sample.
#' @return The function returns the F statistic, Pillai and P value.
#' 
epi_manova <-  function(mixture, model, case_id){
	mixture <- t(mixture)
	mod <- stats::manova(mixture ~ model)
	mod_summary <- summary(mod, tol = 0)$stats
	statistics <- mod_summary[1, c("approx F", "Pillai","Pr(>F)")]
	
	# Calculate the beta mean difference
	keep_case <- which(case_id %in% rownames(mixture))
	case <- mixture[keep_case,]
	controls <- mixture[-keep_case, ]
	coltrols_mean <- colMeans(controls)
	beta_mean_difference <- coltrols_mean - case
	
	output <- list(statistics, beta_mean_difference)
	return(output)
}

#' @title  Creates a data frame containing the results 
#' obtained from MANOVA
#' @description Creates a data frame containing the
#' genomic regions, statistics and direction for the DMRs.
#' @param bump a DMR obtained from \link[bumphunter]{bumphunter}
#' (i.e. a row from \link[bumphunter]{bumphunter} method result).
#' @param beta_bump a beta values matrix for the CpGs in the selected
#' DMR. This matrix is the result of \link[epimutacions]{betas_from_bump}. 
#' @param sts  F statistic, Pillai and P value from
#'  \link[epimutacions]{epi_manova} function results. 
#' @param case a character string specifying the case sample name. 
#' @returns The function returns a data frame containing the following information
#' for each DMR: 
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

res_manova <- function(bump, beta_bump, sts, case) {
	bump$outlier_score <- paste0(sts[[1]][1], "/", sts[[1]][2])
	bump$outlier_significance <- sts[[1]][3]
	bump$adj_pvalue <- NA
	bump$outlier_direction <- ifelse(bump$value < 0, "hypomethylation", "hypermethylation")
	bump$CpG_ids <- paste(rownames(beta_bump), collapse = ",", sep = "")
	bump$sample <- case
	bump[ , c("chr", "start", "end", "sz", "L", "CpG_ids", "outlier_score", "outlier_significance", "adj_pvalue", "outlier_direction", "sample")]
}

filter_manova <- function(bump_out, pvalue_cutoff){
  bump_out$adj_pvalue <- stats::p.adjust(bump_out$outlier_significance, method = "hochberg")   
  bump_out <- bump_out[which(bump_out$adj_pvalue < pvalue_cutoff), , drop = FALSE]
  return(bump_out)
}