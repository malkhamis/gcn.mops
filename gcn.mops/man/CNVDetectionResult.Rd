\name{CNVDetectionResult-class}
\docType{class}

\alias{CNVDetectionResult-class}
\alias{CNVDetectionResult}
\alias{cnvdetectionresult}
\alias{CnvDetectionResult}
\alias{CNVDetectionResult-class}
%\alias{gr,CNVDetectionResult-method}
%\alias{cnvr,CNVDetectionResult-method}
%\alias{cnvs,CNVDetectionResult-method}
%\alias{individualCall,CNVDetectionResult-method}
%\alias{iniCall,CNVDetectionResult-method}
%\alias{integerCopyNumber,CNVDetectionResult-method}
%\alias{localAssessments,CNVDetectionResult-method}
%\alias{normalizedData,CNVDetectionResult-method}
%\alias{params,CNVDetectionResult-method}
%\alias{posteriorProbs,CNVDetectionResult-method}
%\alias{segmentation,CNVDetectionResult-method}
%\alias{segplot,CNVDetectionResult-method}
%\alias{sampleNames,CNVDetectionResult-method}
%\alias{calcIntegerCopyNumbers,CNVDetectionResult-method}
%\alias{calcFractionalCopyNumbers,CNVDetectionResult-method}
%\alias{makeRobustCNVR,CNVDetectionResult-method}


\title{Class "CNVDetectionResult"}
\description{S4 class for storing results of a CNV detection method.}
\section{Slots}{
The following slots are defined for \link{CNVDetectionResult} objects:
\describe{
  \item{gr}{The segments in which the reads are counted. GRanges object.}
  \item{normalizedData}{The normalized data.}
  \item{localAssessments}{The data to which the
  segmentation algorithm is applied. These can be z-Scores,
  ratios, log-ratios or I/NI calls.}
  \item{individualCall}{The CNV call that the method provides for a specific 
  sample}
  \item{iniCall}{The CNV call that the method provides a specific segment.}
  \item{posteriorProbs}{The posterior probabilities for different copy numbers.}
  \item{cnvs}{The detected CNVs.}
  \item{cnvr}{The detected CNV regions.}
  \item{segmentation}{The segmentation of the reference sequence (sample-wise).}
  \item{integerCopyNumber}{The most probable integer copy number.}
  \item{sampleNames}{The sample names.}
  \item{params}{The parameters with which the method was run.}
  }
}  
\section{Methods}{
  \describe{
    \item{gr}{\code{signature(object = "CNVDetectionResult")}: ... }
    \item{cnvr}{\code{signature(object = "CNVDetectionResult")}: ... }
    \item{cnvs}{\code{signature(object = "CNVDetectionResult")}: ... }
    \item{individualCall}{\code{signature(object = "CNVDetectionResult")}: ... }
    \item{iniCall}{\code{signature(object = "CNVDetectionResult")}: ... }
    \item{integerCopyNumber}{\code{signature(object = "CNVDetectionResult")}: ... }
    \item{localAssessments}{\code{signature(object = "CNVDetectionResult")}: ... }
    \item{normalizedData}{\code{signature(object = "CNVDetectionResult")}: ... }
    \item{params}{\code{signature(object = "CNVDetectionResult")}: ... }
    \item{plot}{\code{signature(x = "CNVDetectionResult", y = "missing")}: ... }
    \item{posteriorProbs}{\code{signature(object = "CNVDetectionResult")}: ... }
    \item{sampleNames}{\code{signature(object = "CNVDetectionResult")}: ... }
    \item{segmentation}{\code{signature(object = "CNVDetectionResult")}: ... }
    \item{segplot}{\code{signature(object = "CNVDetectionResult")}: ... }
    \item{show}{\code{signature(object = "CNVDetectionResult")}: ... }
    \item{calcIntegerCopyNumbers}{\code{signature(object = "CNVDetectionResult")}: ... }
    \item{makeRobustCNVR}{\code{signature(object = "CNVDetectionResult")}: ... }
	 }
 }
\author{
  Guenter Klambauer \email{klambauer@bioinf.jku.at}
}

\examples{
showClass("CNVDetectionResult")
}
\keyword{classes}