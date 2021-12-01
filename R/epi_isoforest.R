#' @title Identifies epimutations using Isolation Forest
#' @description  This function identifies regions with CpGs being outliers
#' using \link[isotree]{isolation.forest} approach. 
#' @param mixture beta values matrix. Samples in columns and
#' CpGs in rows.
#' @param case_id a character string specifying the name of the case sample.
#' @param ntrees number of binary trees to build for the model. 
#' Default is \code{100}. 
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
	score <- stats::predict(iso, test)
	
	
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
#' @param outlier_score_cutoff numeric specifying the outlier score cut off
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

res_isoforest <- function(bump, sts, outlier_score_cutoff){
  if(sts > outlier_score_cutoff){
	bump$outlier_score <- sts
	bump$pvalue <- NA
	bump$adj_pvalue <- NA
	bump$outlier_direction <- ifelse(bump$value < 0, "hypomethylation", 
	                                 "hypermethylation")
	bump[ , c("chromosome", 
	          "start", "end",
	          "sz", "cpg_n", 
	          "cpg_ids", "outlier_score",
	          "outlier_direction", "pvalue", 
	          "adj_pvalue",  "sample")]
	}else{
    data.frame(chromosome = character(), start = numeric(), 
               end = numeric(),
               sz = numeric(),
               cpg_n = numeric(),
               cpg_ids = character(),
               outlier_score = numeric(),
               outlier_direction = character(),
               pvalue = numeric(),
               adj_pvalue = numeric(),
               sample = character())
  }
	 }