\name{segmentation}
\alias{segmentation}
\title{This generic function returns segmentation of
 a CNV detection method stored in an instance of
\code{\link{CNVDetectionResult-class}}.}
\arguments{
  \item{object}{An instance of "CNVDetectionResult"}
}
\value{
  \code{segmentation} returns a eturns a "GRanges" object
  containing the segmentation.
}
\description{
  This generic function returns segmentation of a CNV
  detection method stored in an instance of
  \code{\link{CNVDetectionResult-class}}.
}
\examples{
data(cn.mops)
r <- cn.mops(X[1:100,1:5])
segmentation(r)
}
\author{
  Guenter Klambauer \email{klambauer@bioinf.jku.at}
}

