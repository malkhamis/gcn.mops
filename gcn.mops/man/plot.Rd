\name{plot}
\docType{methods}
\alias{plot-methods}
\alias{plot,CNVDetectionResult,missing-method}

\title{Plots a CNVDetectionResult}
\description{
  Plots read counts, call values and CNV calls in an
  identified CNV region.
}
\usage{
\S4method{plot}{CNVDetectionResult,missing}(x,
             which,margin=c(10,10),toFile=FALSE)
}
\arguments{
  \item{x}{An instance of "CNVDetectionResult"}

  \item{which}{The index of the CNV region to be plotted.}

  \item{margin}{Vector of two positive integers that states
  how many segments left and right of the CNV region should
  be included in the plot. Default = c(10,10).}

  \item{toFile}{Logical value whether the output should be
  plotted to a file. Default = FALSE.}
}
\value{
  Generates a CNV calling plot.
}

\author{
  Guenter Klambauer \email{klambauer@bioinf.jku.at}
}

