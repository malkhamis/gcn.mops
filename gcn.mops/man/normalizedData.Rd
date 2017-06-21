\name{normalizedData}
\alias{normalizedData}
\title{This generic function returns the normalized data
 of a CNV detection method stored in an instance of
\code{\link{CNVDetectionResult-class}}.}
\arguments{
  \item{object}{An instance of "CNVDetectionResult".}
}
\value{
  \code{normalizedData} returns a "GRanges" object
  containing the normalized data.
}
\description{
  This generic function returns the normalized data of a
  CNV detection method stored in an instance of
  \code{\link{CNVDetectionResult-class}}.
}
\examples{
data(cn.mops)
r <- cn.mops(X[1:100,1:5])
normalizedData(r)
}
\author{
  Guenter Klambauer \email{klambauer@bioinf.jku.at}
}

