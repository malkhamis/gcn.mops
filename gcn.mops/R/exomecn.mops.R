#' @title Copy number detection in exome sequencing data.
#' @description Performs the cn.mops algorithm for copy number detection in
#' NGS data with parameters adjusted to exome sequencing data.
#' 
#' @param input Either an instance of "GRanges" or a raw data matrix, where
#' columns are interpreted as samples and rows as genomic regions. An entry is
#' the read count of a sample in the genomic region.
#' @param I Vector positive real values that contain the expected fold change
#' of the copy number classes.  Length of this vector must be equal to the 
#' length of the "classes" parameter vector. For human copy number polymorphisms 
#' we suggest to use the default I = c(0.025,0.5,1,1.5,2,2.5,3,3.5,4).
#' @param classes Vector of characters of the same length as the parameter
#' vector "I". One vector element must be named "CN2". The names reflect the 
#' labels of the copy number classes. 
#' Default = c("CN0","CN1","CN2","CN3","CN4","CN5","CN6","CN7","CN8").
#' @param priorImpact Positive real value that reflects how strong the prior
#' assumption affects the result. The higher the value the more samples will
#' be assumed to have copy number 2. Default = 10.
#' @param cyc Positive integer that sets the number of cycles for the algorithm.
#' Usually after less than 15 cycles convergence is reached. Default = 20.
#' @param parallel How many cores are used for the computation. If set to zero
#' than no parallelization is applied. Default = 0.
#' @param norm The normalization strategy to be used. 
#' If set to 0 the read counts are not normalized and cn.mops does not model 
#' different coverages. 
#' If set to 1 the read counts are normalized. 
#' If set to 2 the read counts are not normalized and cn.mops models different
#' coverages. (Default=1).
#' @param normType Mode of the normalization technique. Possible values are 
#' "mean","min","median","quant", "poisson" and "mode". 
#' Read counts will be scaled sample-wise. Default = "poisson".
#' @param sizeFactor  By this parameter one can decide to how the size factors 
#' are calculated.
#' Possible choices are the the mean, median or mode coverage ("mean", "median", "mode") or any quantile 
#' ("quant").
#' @param normQu Real value between 0 and 1.  
#' If the "normType" parameter is set to "quant" then this parameter sets the 
#' quantile that is used for the normalization. Default = 0.25. 
#' @param quSizeFactor Quantile of the sizeFactor if sizeFactor is set to "quant".
#' 0.75 corresponds to "upper quartile normalization". Real value between 0 and 1. Default = 0.75.
#' @param upperThreshold Positive real value that sets the cut-off for copy
#' number gains. All CNV calling values above this value will be called as 
#' "gain". The value should be set close to the log2 of the expected foldchange
#' for copy number 3 or 4. Default = 0.55.
#' @param lowerThreshold Negative real value that sets the cut-off for copy
#' number losses. All CNV calling values below this value will be called as 
#' "loss". The value should be set close to the log2 of the expected foldchange
#' for copy number 1 or 0. Default = -0.8.
#' @param minWidth Positive integer that is exactly the parameter "min.width"
#' of the "segment" function of "DNAcopy". minWidth is the minimum number 
#' of segments a CNV should span. Default = 5.
#' @param segAlgorithm Which segmentation algorithm should be used. If set to
#' "DNAcopy" circular binary segmentation is performed. Any other value will
#' initiate the use of our fast segmentation algorithm. Default = "fast".
#' @param minReadCount If all samples are below this value the algorithm will
#' return the prior knowledge. This prevents that the algorithm from being 
#' applied to segments with very low coverage. Default=1. 
#' @param useMedian Whether "median" instead of "mean" of a segment
#' should be used for the CNV call. Default=FALSE. 
#' @param returnPosterior Flag that decides whether the posterior probabilities
#' should be returned. The posterior probabilities have a dimension of samples
#' times copy number states times genomic regions and therefore consume a lot
#' of memory. Default=FALSE.
#' @param ... Additional parameters will be passed to the "DNAcopy"
#' or the standard segmentation algorithm.
#' @examples 
#' data(cn.mops)
#' exomecn.mops(exomeCounts)
#'
#' @useDynLib cn.mops
#' @return An instance of "CNVDetectionResult".
#' @author Guenter Klambauer \email{klambauer@@bioinf.jku.at}
#' @export


exomecn.mops <- function(input,I = c(0.025,0.5,1,1.5,2,2.5,3,3.5,4),
		classes=c("CN0","CN1","CN2","CN3","CN4","CN5","CN6","CN7","CN8"),
		priorImpact = 10,cyc = 20,parallel=0,
		norm=1, normType="poisson",sizeFactor="mean",normQu=0.25, quSizeFactor=0.75,
		upperThreshold=0.5,lowerThreshold=-0.8,
		minWidth=5,segAlgorithm="fast",minReadCount=1,useMedian=FALSE,returnPosterior=FALSE,...){
	res <- cn.mops(input=input,I=I,classes=classes,priorImpact=priorImpact,
			cyc=cyc,parallel=parallel,normType=normType,normQu=normQu,norm=norm,sizeFactor=sizeFactor,
			quSizeFactor=quSizeFactor, lowerThreshold=lowerThreshold,
			upperThreshold=upperThreshold,minWidth=minWidth,segAlgorithm=segAlgorithm,
			minReadCount=minReadCount,useMedian=useMedian,returnPosterior=returnPosterior,...)
	return(res)
}
	