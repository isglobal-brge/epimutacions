#' @title Identifies epimutations using Isolation Forest
#' @description  This function identifies regions with CpGs being outliers
#' using \link[isotree]{isolation.forest} approach. 
#' @param mixture beta values matrix. Samples in columns and
#' CpGs in rows.
#' @param case_id a character string specifying the name of the case sample.
#' @return The function returns the outlier score for the given case sample.
#' 
epi_isoforest <- function(mixture, case_id, ntrees) {
	mixture <- t(mixture)
	#Generate train and test(sample with suspected disease) data frame
	train <- mixture[row.names(mixture) != case_id, ]
	test <- mixture[case_id, , drop = FALSE]
	
	#Run the isolation forest methods 
	iso <- isotree::isolation.forest(train, ntrees = ntrees)
	
	#Predict
	score <- predict(iso, test)
	
	return(score)
}

#' @title  Creates a data frame containing the results 
#' obtained from Isolation Forest
#' @description Creates a data frame containing the
#' genomic regions, statistics and direction for the DMRs.
#' @param bump a DMR obtained from \link[bumphunter]{bumphunter}
#' (i.e. a row from \link[bumphunter]{bumphunter} method result).
#' @param beta_bump a beta values matrix for the CpGs in the selected
#' DMR. This matrix is the result of \link[epimutacions]{betas_from_bump}. 
#' @param sts the outlier score from
#'  \link[epimutacions]{epi_isoforest} function results. 
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

res_isoforest <- function(bump, beta_bump, sts, case) {
	bump$outlier_score <- sts
	bump$outlier_significance <- NA
	bump$outlier_direction <- ifelse(bump$value < 0, "hypomethylation", "hypermethylation")
	bump$CpG_ids <- paste(rownames(beta_bump), collapse = ",", sep = "")
	bump$sample <- case
	bump[ , c("chr", "start", "end", "sz", "L", "CpG_ids", "outlier_score", "outlier_significance", "outlier_direction", "sample")]
}