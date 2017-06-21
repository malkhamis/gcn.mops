# Copyright (C) 2011 Klambauer Guenter 
# <klambauer@bioinf.jku.at>
#'  This generic function returns the genomic ranges
#'  of a CNV detection method stored in an instance of 
#' \code{\link{CNVDetectionResult-class}}.
#' 
#' @param object An instance of "CNVDetectionResult".
#' @examples
#' data(cn.mops)
#' r <- cn.mops(X[1:100,1:5])
#' gr(r)
#' @return \code{normalizedData} returns a "GRanges" object containing
#' the normalized data.
#' @author Guenter Klambauer \email{klambauer@@bioinf.jku.at}
#' @export

setMethod("gr", signature = "CNVDetectionResult", 
		definition = function(object) object@gr)



#'  This generic function returns the normalized data
#'  of a CNV detection method stored in an instance of 
#' \code{\link{CNVDetectionResult-class}}.
#' 
#' @param object An instance of "CNVDetectionResult".
#' @examples
#' data(cn.mops)
#' r <- cn.mops(X[1:100,1:5])
#' normalizedData(r)
#' @return \code{normalizedData} returns a "GRanges" object containing
#' the normalized data.
#' @author Guenter Klambauer \email{klambauer@@bioinf.jku.at}
#' @export

setMethod("normalizedData", signature = "CNVDetectionResult", 
		definition = function(object){
			xx <- object@gr
			values(xx) <- object@normalizedData
			return(xx)
		}
)



#'  This generic function returns the local assessments, i.e. 
#' signed individual informative/non-informative calls,
#' of a CNV detection method stored in an instance of 
#' \code{\link{CNVDetectionResult-class}}. For other CNV detection methods
#' this can be (log-) ratios or z-scores.
#' 
#' @param object An instance of "CNVDetectionResult" 
#' @examples
#' data(cn.mops)
#' r <- cn.mops(X[1:100,1:5])
#' localAssessments(r)
#' @return \code{localAssessments} returns a "GRanges" object containing
#' the local assessments.
#' @author Guenter Klambauer \email{klambauer@@bioinf.jku.at}
#' @export

setMethod("localAssessments", signature = "CNVDetectionResult", 
		definition = function(object){
			xx <- object@gr
			values(xx) <- object@localAssessments
			return(xx)
		}
)




#'  This generic function returns the individual calls of
#'  a CNV detection method stored in an instance of 
#' \code{\link{CNVDetectionResult-class}}. 
#' 
#' @param object An instance of "CNVDetectionResult" 
#' @examples
#' data(cn.mops)
#' r <- cn.mops(X[1:100,1:5])
#' individualCall(r)
#' @return \code{individualCalls} returns a "GRanges" object containing
#' the individual calls.
#' @author Guenter Klambauer \email{klambauer@@bioinf.jku.at}
#' @export

setMethod("individualCall", signature = "CNVDetectionResult", 
		definition = function(object){
			xx <- object@gr
			values(xx) <- object@individualCall
			return(xx)
		}
)




#'  This generic function returns the informative/non-informative call of
#'  a CNV detection method stored in an instance of 
#' \code{\link{CNVDetectionResult-class}}.
#'  The I/NI call is a measure for a genomic
#' segment across all samples, whether this segment is a CNV region 
#' (informative) or a normal genomic region (non-informative).
#' 
#' @param object An instance of "CNVDetectionResult" 
#' @examples
#' data(cn.mops)
#' r <- cn.mops(X[1:100,1:5])
#' iniCall(r)
#' @return \code{iniCall} returns a "GRanges" object containing
#' the individual calls.
#' @author Guenter Klambauer \email{klambauer@@bioinf.jku.at}
#' @export


setMethod("iniCall", signature = "CNVDetectionResult", 
		definition = function(object){
			xx <- object@gr
			values(xx) <- object@iniCall
			return(xx)
		}
)


#'  This generic function returns the posterior probabilities of
#'  a CNV detection method stored in an instance of 
#' \code{\link{CNVDetectionResult-class}}. 
#' The posterior probabilities are represented
#' as a three dimensional array, where the three dimensions are segment,
#' copy number and individual.
#' 
#' @param object An instance of "CNVDetectionResult" 
#' @examples
#' data(cn.mops)
#' r <- cn.mops(X[1:100,1:5])
#' posteriorProbs(r)
#' @return \code{posteriorProbs} returns a  three dimensional array.
#' @author Guenter Klambauer \email{klambauer@@bioinf.jku.at}
#' @export

setMethod("posteriorProbs", signature = "CNVDetectionResult", 
		definition = function(object) object@posteriorProbs)



#'  This generic function returns CNVs of
#'  a CNV detection method stored in an instance of 
#' \code{\link{CNVDetectionResult-class}}. 
#' 
#' @param object An instance of "CNVDetectionResult" 
#' @examples
#' data(cn.mops)
#' r <- cn.mops(X[1:100,1:5])
#' cnvs(r)
#' @return \code{cnvs} returns a  eturns a "GRanges" object containing
#' the CNVs.
#' @author Guenter Klambauer \email{klambauer@@bioinf.jku.at}
#' @export

setMethod("cnvs", signature = "CNVDetectionResult", 
		definition = function(object) object@cnvs)




#'  This generic function returns CNV regions of
#'  a CNV detection method stored in an instance of 
#' \code{\link{CNVDetectionResult-class}}. 
#' 
#' @param object An instance of "CNVDetectionResult" 
#' @examples
#' data(cn.mops)
#' r <- cn.mops(X[1:100,1:5])
#' cnvr(r)
#' @return \code{cnvr} returns a  eturns a "GRanges" object containing
#' the CNV regions.
#' @author Guenter Klambauer \email{klambauer@@bioinf.jku.at}
#' @export

setMethod("cnvr", signature = "CNVDetectionResult", 
		definition = function(object) object@cnvr)




#'  This generic function returns segmentation of
#'  a CNV detection method stored in an instance of 
#' \code{\link{CNVDetectionResult-class}}. 
#' 
#' @param object An instance of "CNVDetectionResult" 
#' @examples
#' data(cn.mops)
#' r <- cn.mops(X[1:100,1:5])
#' segmentation(r)
#' @return \code{segmentation} returns a  eturns a "GRanges" object containing
#' the segmentation.
#' @author Guenter Klambauer \email{klambauer@@bioinf.jku.at}
#' @export

setMethod("segmentation", signature = "CNVDetectionResult", 
		definition = function(object) object@segmentation)




#'  This generic function returns the integer copy numbers of
#'  a CNV detection method stored in an instance of 
#' \code{\link{CNVDetectionResult-class}}. 
#' 
#' @param object An instance of "CNVDetectionResult" 
#' @examples
#' data(cn.mops)
#' r <- cn.mops(X[1:100,1:5])
#' integerCopyNumber(r)
#' @return \code{integerCopyNumber} returns a  eturns a "GRanges" object containing
#' the integer copy numbers.
#' @author Guenter Klambauer \email{klambauer@@bioinf.jku.at}
#' @export



setMethod("integerCopyNumber", signature = "CNVDetectionResult", 
		definition = function(object){
			xx <- object@gr
			values(xx) <- object@integerCopyNumber
			return(xx)
		}
)



#'  This generic function returns the parameters of
#'  a CNV detection method stored in an instance of 
#' \code{\link{CNVDetectionResult-class}}. 
#' 
#' @param object An instance of "CNVDetectionResult" 
#' @examples
#' data(cn.mops)
#' r <- cn.mops(X[1:100,1:5])
#' params(r)
#' @return \code{params} returns a  eturns a "GRanges" object containing
#' the parameters.
#' @author Guenter Klambauer \email{klambauer@@bioinf.jku.at}
#' @export

setMethod("params", signature = "CNVDetectionResult", 
		definition = function(object) object@params)




#'  This generic function returns the sample names of
#'  a CNV detection method stored in an instance of 
#' \code{\link{CNVDetectionResult-class}}. 
#' 
#' @param object An instance of "CNVDetectionResult" 
#' @examples
#' data(cn.mops)
#' r <- cn.mops(X[1:100,1:5])
#' sampleNames(r)
#' @return \code{sampleNames} returns a  eturns a "GRanges" object containing
#' the parameters.
#' @author Guenter Klambauer \email{klambauer@@bioinf.jku.at}
#' @export
setMethod("sampleNames", signature = "CNVDetectionResult", 
		definition = function(object) object@sampleNames)


















