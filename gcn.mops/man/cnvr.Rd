\name{cnvr}
\alias{cnvr}
\title{This generic function returns CNV regions of
 a CNV detection method stored in an instance of
\code{\link{CNVDetectionResult-class}}.}
\arguments{
  \item{object}{An instance of "CNVDetectionResult"}
}
\value{
  \code{cnvr} returns a eturns a "GRanges" object
  containing the CNV regions.
}
\description{
  This generic function returns CNV regions of a CNV
  detection method stored in an instance of
  \code{\link{CNVDetectionResult-class}}.
}
\examples{
data(cn.mops)
r <- cn.mops(X[1:100,1:5])
cnvr(r)
}
\author{
  Guenter Klambauer \email{klambauer@bioinf.jku.at}
}

