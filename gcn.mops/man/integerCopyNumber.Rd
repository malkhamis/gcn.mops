\name{integerCopyNumber}
\alias{integerCopyNumber}
\title{This generic function returns the integer copy numbers of
 a CNV detection method stored in an instance of
\code{\link{CNVDetectionResult-class}}.}
\arguments{
  \item{object}{An instance of "CNVDetectionResult"}
}
\value{
  \code{integerCopyNumber} returns a eturns a "GRanges"
  object containing the integer copy numbers.
}
\description{
  This generic function returns the integer copy numbers of
  a CNV detection method stored in an instance of
  \code{\link{CNVDetectionResult-class}}.
}
\examples{
data(cn.mops)
r <- cn.mops(X[1:100,1:5])
integerCopyNumber(r)
}
\author{
  Guenter Klambauer \email{klambauer@bioinf.jku.at}
}

