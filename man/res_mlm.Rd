% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/epi_mlm.R
\name{res_mlm}
\alias{res_mlm}
\title{Creates a data frame containing the results
obtained from MLM}
\usage{
res_mlm(bump, sts)
}
\arguments{
\item{bump}{a DMR obtained from \link[bumphunter]{bumphunter}
(i.e. a row from \link[bumphunter]{bumphunter} method result).}

\item{sts}{the F statistic, R2 test statistic
and Pillai obtained as a result
of \link[epimutacions]{epi_mlm} function.}
}
\value{
The function returns a data
frame containing the following
information for each DMR:
\itemize{
\item genomic ranges
\item DMR base pairs
\item number and name of CpGs in DMR
\item statistics:
\itemize{
\item Outlier score
\item Outlier significance
\item Outlier direction
}
\item Sample name
}

For more information about the output see
\link[epimutacions]{epimutations}.
}
\description{
Creates a data frame containing the
genomic regions, statistics and direction for the DMRs.
}
