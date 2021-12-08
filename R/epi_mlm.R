#' @title  Detects epimutations using Multivariate Linear Model (MLM)
#' @description  Identifies CpGs with outlier methylation values 
#' using methylated Multivariate Linear Model
#' @param mixture beta values matrix.  Samples in columns and
#' CpGs in rows. 
#' @param  model design (or model) matrix.
#' @return The function returns the F statistic, R2 test statistic and Pillai.
#' 
epi_mlm <- function(mixture, model) {
  
  mod <- mlm(t(mixture) ~ model[,2])
	statistics <- mod$aov.tab[1, c("F value", "R2", "Pr(>F)")]
	return(statistics)
}

#' @title  Creates a data frame containing the results 
#' obtained from MLM
#' @description Creates a data frame containing the
#' genomic regions, statistics and direction for the DMRs.
#' @param bump a DMR obtained from \link[bumphunter]{bumphunter}
#' (i.e. a row from \link[bumphunter]{bumphunter} method result).
#' @param sts the F statistic, R2 test statistic and Pillai obtained as a result
#' of \link[epimutacions]{epi_mlm} function. 
#' @returns The function returns a data frame containing the following 
#' information for each DMR: 
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
#' 

res_mlm <- function(bump, sts) {
	bump$outlier_score <- paste0(sts[1], "/", sts[2])
	bump$outlier_direction <- ifelse(bump$value < 0, "hypomethylation", 
	                                 "hypermethylation")
	bump$pvalue <- sts[3]
	bump$adj_pvalue <- NA
	bump[ , c("chromosome", "start", "end", "sz", "cpg_n", 
	          "cpg_ids", "outlier_score",
	          "outlier_direction", "pvalue", 
	          "adj_pvalue",  "sample")]
}
