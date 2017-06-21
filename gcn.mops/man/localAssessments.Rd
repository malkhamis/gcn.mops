\name{localAssessments}
\alias{localAssessments}
\title{This generic function returns the local assessments, i.e.
signed individual informative/non-informative calls,
of a CNV detection method stored in an instance of
\code{\link{CNVDetectionResult-class}}. For other CNV detection methods
this can be (log-) ratios or z-scores.}
\arguments{
  \item{object}{An instance of "CNVDetectionResult"}
}
\value{
  \code{localAssessments} returns a "GRanges" object
  containing the local assessments.
}
\description{
  This generic function returns the local assessments, i.e.
  signed individual informative/non-informative calls, of a
  CNV detection method stored in an instance of
  \code{\link{CNVDetectionResult-class}}. For other CNV
  detection methods this can be (log-) ratios or z-scores.
}
\examples{
data(cn.mops)
r <- cn.mops(X[1:100,1:5])
localAssessments(r)
}
\author{
  Guenter Klambauer \email{klambauer@bioinf.jku.at}
}

