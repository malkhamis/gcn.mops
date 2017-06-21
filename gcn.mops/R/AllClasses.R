# Copyright (C) 2011 Klambauer Guenter 
# <klambauer@bioinf.jku.at>

# S4 class definition for the result object of a CNV detection method

#' @export

setClass("CNVDetectionResult",
		representation = representation
				(
				gr					= "GRanges",
				normalizedData     	= "matrix",
				localAssessments    = "matrix",
				individualCall      = "matrix",
				iniCall        		= "numeric",
				posteriorProbs		= "array",
				cnvs				= "GRanges",
				cnvr				= "GRanges",
				segmentation		= "GRanges",
				integerCopyNumber	= "matrix",
				params				= "list",
				sampleNames			= "character"
		),
		prototype = prototype
				(
				gr					= GRanges(),
				normalizedData     	= matrix(),
				localAssessments    = matrix(),
				individualCall      = matrix(),
				iniCall        		= vector("numeric",1),
				posteriorProbs		= array(NA,dim=c(1,1,1)),
				cnvs				= GRanges(),
				cnvr				= GRanges(),
				segmentation		= GRanges(),
				integerCopyNumber	= matrix(),
				params				= list(),
				sampleNames			= vector("character",1)
				)
)





