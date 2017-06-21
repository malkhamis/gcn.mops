\name{sampleNames}
\alias{sampleNames}
\title{This generic function returns the sample names of
 a CNV detection method stored in an instance of
\code{\link{CNVDetectionResult-class}}.}
\arguments{
  \item{object}{An instance of "CNVDetectionResult"}
}
\value{
  \code{sampleNames} returns a eturns a "GRanges" object
  containing the parameters.
}
\description{
  This generic function returns the sample names of a CNV
  detection method stored in an instance of
  \code{\link{CNVDetectionResult-class}}.
}
\examples{
data(cn.mops)
r <- cn.mops(X[1:100,1:5])
sampleNames(r)
}
\author{
  Guenter Klambauer \email{klambauer@bioinf.jku.at}
}

