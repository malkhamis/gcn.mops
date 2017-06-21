\name{calcIntegerCopyNumbers}
\alias{calcIntegerCopyNumbers}
\title{Calculation of integer copy numbers for the CNVs and CNV regions.}
\arguments{
  \item{object}{An instance of "CNVDetectionResult"}
}
\value{
  \code{calcIntegerCopyNumbers} returns an instance of
  "CNVDetectionResult".
}
\description{
  This generic function calculates the integer copy numbers
  of a CNV detection method stored in an instance of
  \code{\link{CNVDetectionResult-class}}.
}
\examples{
data(cn.mops)
r <- cn.mops(X[1:100,1:5])
calcIntegerCopyNumbers(r)
}
\author{
  Guenter Klambauer \email{klambauer@bioinf.jku.at}
}

