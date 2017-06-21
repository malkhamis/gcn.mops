\name{params}
\alias{params}
\title{This generic function returns the parameters of
 a CNV detection method stored in an instance of
\code{\link{CNVDetectionResult-class}}.}
\arguments{
  \item{object}{An instance of "CNVDetectionResult"}
}
\value{
  \code{params} returns a eturns a "GRanges" object
  containing the parameters.
}
\description{
  This generic function returns the parameters of a CNV
  detection method stored in an instance of
  \code{\link{CNVDetectionResult-class}}.
}
\examples{
data(cn.mops)
r <- cn.mops(X[1:100,1:5])
params(r)
}
\author{
  Guenter Klambauer \email{klambauer@bioinf.jku.at}
}

