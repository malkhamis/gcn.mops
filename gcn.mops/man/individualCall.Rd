\name{individualCall}
\alias{individualCall}
\title{This generic function returns the individual calls of
 a CNV detection method stored in an instance of
\code{\link{CNVDetectionResult-class}}.}
\arguments{
  \item{object}{An instance of "CNVDetectionResult"}
}
\value{
  \code{individualCalls} returns a "GRanges" object
  containing the individual calls.
}
\description{
  This generic function returns the individual calls of a
  CNV detection method stored in an instance of
  \code{\link{CNVDetectionResult-class}}.
}
\examples{
data(cn.mops)
r <- cn.mops(X[1:100,1:5])
individualCall(r)
}
\author{
  Guenter Klambauer \email{klambauer@bioinf.jku.at}
}

