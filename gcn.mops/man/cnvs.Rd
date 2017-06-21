\name{cnvs}
\alias{cnvs}
\title{This generic function returns CNVs of
 a CNV detection method stored in an instance of
\code{\link{CNVDetectionResult-class}}.}
\arguments{
  \item{object}{An instance of "CNVDetectionResult"}
}
\value{
  \code{cnvs} returns a eturns a "GRanges" object
  containing the CNVs.
}
\description{
  This generic function returns CNVs of a CNV detection
  method stored in an instance of
  \code{\link{CNVDetectionResult-class}}.
}
\examples{
data(cn.mops)
r <- cn.mops(X[1:100,1:5])
cnvs(r)
}
\author{
  Guenter Klambauer \email{klambauer@bioinf.jku.at}
}

