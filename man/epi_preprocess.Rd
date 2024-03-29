% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/epi_preprocess.R
\name{epi_preprocess}
\alias{epi_preprocess}
\title{Preprocess methylation array}
\usage{
epi_preprocess(
  cases_dir,
  reference_panel,
  pattern = "csv$",
  normalize = "raw",
  norm_param = norm_parameters(),
  verbose = FALSE
)
}
\arguments{
\item{cases_dir}{the base directory from which the search is started.}

\item{reference_panel}{an \link[minfi]{RGChannelSet}
object containing the reference
panel (controls) samples.}

\item{pattern}{What pattern is used to identify a sample sheet file.}

\item{normalize}{a character string
specifying the selected preprocess method.
For more information see \strong{Details} or
\href{shorturl.at/nwIN2}{minfi package user's Guide}.
It can be set as: \code{"raw"}, \code{"illumina"},
\code{"swan"}, \code{"quantile"},
\code{"noob"} or \code{"funnorm"}.)}

\item{norm_param}{the parameters for each preprocessing method.
See the function \link[epimutacions]{norm_parameters}.}

\item{verbose}{logical. If TRUE additional
details about the procedure will provide to the user.
The default is FALSE.}
}
\value{
\code{epi_preprocess} function returns a
\link[minfi]{GenomicRatioSet} object
containing case and control (reference panel) samples.
}
\description{
The \code{epi_preprocess} function reads
Illumina methylation sample
sheet for case samples and it merges them with
\link[minfi]{RGChannelSet} reference panel.
The final dataset is normalized using minfi package preprocess methods.
}
\details{
The \code{epi_preprocess} function reads Illumina methylation sample
sheet for case samples and it merges them with
\link[minfi]{RGChannelSet} reference panel.
The final dataset is normalized using
different minfi package preprocess methods:
\itemize{
\item \code{"raw"}: \link[minfi]{preprocessRaw}
\item \code{"illumina"}: \link[minfi]{preprocessIllumina}
\item \code{"swan"}: \link[minfi]{preprocessSWAN}
\item \code{"quantile"}: \link[minfi]{preprocessQuantile}
\item \code{"noob"}: \link[minfi]{preprocessNoob}
\item \code{"funnorm"}: \link[minfi]{preprocessFunnorm}
}
}
\examples{

# The reference panel for this example is available in 
#epimutacionsData (ExperimentHub) package
\donttest{
   library(ExperimentHub)
   eh <- ExperimentHub()
   query(eh, c("epimutacionsData"))
   reference_panel <- eh[["EH6691"]]
   cases_dir <- system.file("extdata", package = "epimutacionsData")
   #Preprocessing

   epi_preprocess( cases_dir, 
                   reference_panel, 
                   pattern = "SampleSheet.csv")
}


}
