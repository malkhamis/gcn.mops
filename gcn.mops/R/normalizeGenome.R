# Copyright (C) 2011 Klambauer Guenter 
# <klambauer@bioinf.jku.at>


#' @title Normalization of NGS data
#' 
#' @description Normalize quantitative NGS data in order to make counts comparable over
#' samples. Scales each samples' reads such that the coverage is even for
#' all samples after normalization. 
#' @param X Matrix of positive real values, where
#' columns are interpreted as samples and rows as genomic regions. An entry is
#' the read count of a sample in the genomic region.  Alternatively this can be
#' a GRanges object containing the read counts as values.
#' @param normType Type of the normalization technique. Each samples'
#' read counts are scaled such that the total number of reads are comparable across 
#' samples.
#' If this parameter is set to the value "mode", 
#' the read counts are scaled such that each samples'
#' most frequent value (the "mode") is equal after normalization. 
#' Accordingly for the other options are "mean","median","poisson", "quant", and "mode". 
#' Default = "poisson".
#' @param sizeFactor  By this parameter one can decide to how the size factors 
#' are calculated.
#' Possible choices are the the mean, median or mode coverage ("mean", "median", "mode") or any quantile 
#' ("quant").
#' @param qu Quantile of the normType if normType is set to "quant" .Real value between 0 and 1. Default = 0.25.
#' @param quSizeFactor Quantile of the sizeFactor if sizeFactor is set to "quant".
#' 0.75 corresponds to "upper quartile normalization". Real value between 0 and 1. Default = 0.75.
#' @param ploidy An integer value for each sample or each column in the read
#' count matrix. At least two samples must have a ploidy of 2. Default = "missing".
#' @examples 
#' data(cn.mops)
#' X.norm <- normalizeGenome(X)
#' @return A data matrix of normalized read counts with the same dimensions 
#' as the input matrix X.
#' @author Guenter Klambauer \email{klambauer@@bioinf.jku.at}
#' @export



normalizeGenome <- function(X,normType="poisson", sizeFactor="mean", qu=0.25, quSizeFactor=0.75, ploidy){
	if (class(X)=="GRanges"){
		X.counts <- as.matrix(values(X))
		chr <- rep("Chr",nrow(X.counts))
		YY <- normalizeChromosomes(X.counts, chr=chr, normType=normType, 
				sizeFactor=sizeFactor, qu=qu, quSizeFactor=quSizeFactor, ploidy=ploidy)
		values(X) <- YY
		return(X)
	} else {
		YY <- normalizeChromosomes(X,normType=normType, 
				sizeFactor=sizeFactor, qu=qu, quSizeFactor=quSizeFactor, ploidy=ploidy)
		return(YY)
	}
}
