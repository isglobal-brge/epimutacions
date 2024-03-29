% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/betas_from_bump.R
\name{betas_from_bump}
\alias{betas_from_bump}
\title{Obtains bumps beta values}
\usage{
betas_from_bump(bump, fd, betas)
}
\arguments{
\item{bump}{the result from \link[bumphunter]{bumphunter}.}

\item{fd}{a data frame containing the genomic ranges for each CpGs.}

\item{betas}{a matrix containing the beta
values for all CpGs in each sample.}
}
\value{
The function returns a data
frame containing the beta values
for each sample and CpG into DMR.
}
\description{
The function obtains beta values
corresponding to the CpGs into DMRs.
}
