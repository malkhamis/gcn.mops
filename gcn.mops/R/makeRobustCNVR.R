#' This generic function calculates robust CNV regions by 
#' segmenting the I/NI call per genomic region
#' of an object \code{\link{CNVDetectionResult-class}}.
#' 
#' @name makeRobustCNVR
#' @title Calculates robust CNV regions.
#' @description This generic function calculates robust CNV regions by 
#' segmenting the I/NI call per genomic region
#' of an object \code{\link{CNVDetectionResult-class}}.
#' @aliases makeRobustCNVR,CNVDetectionResult-method
#' @param object An instance of "CNVDetectionResult" 
#' @param robust Robustness parameter. The higher the value, the more samples
#' are required to have a CNV that confirms the CNV region. Setting this 
#' parameter to 0 restores the original CNV regions.
#' (Default=0.5)
#' @param minWidth The minimum length measured in genomic regions a CNV region
#' has to span in order to be called. A parameter of the segmentation algorithm.
#' (Default=4).
#' @param ... Additional parameters passed to the segmentation algorithm.
#' @details cn.mops usually reports a CNV region if at least one individual
#' has a CNV in this region. For some applications it is useful to find more
#' common CNV regions, i.e., regions in which more than one sample has a CNV.
#' The I/NI call measures both signal strength and how many sample show an
#' abnormal copy number, therefore segmentation of the I/NI call can provide
#' robust CNV regions. 
#' 
#' @examples
#' data(cn.mops)
#' r <- cn.mops(X[1:100,1:5])
#' rr <- calcIntegerCopyNumbers(makeRobustCNVR(r,robust=0.1,minWidth=3))
#' @return \code{\link{makeRobustCNVR}} returns a "CNVDetectionResult" 
#' object containing new values in the slot "cnvr".
#' @author Guenter Klambauer \email{klambauer@@bioinf.jku.at}
#' @export



setMethod("makeRobustCNVR", signature = "CNVDetectionResult", 
		definition = function(object, robust=0.5, minWidth=4,...) {
			
			if (robust>0){
				chr <- as.character(seqnames(object@gr))
				chrOrder <- unique(chr) 
				chrBpts <- cumsum(table(chr)[chrOrder])
				chrDf <- data.frame(start=c(1,chrBpts[-length(chrBpts)]+1),
						end=chrBpts)
				rownames(chrDf) <- chrOrder
				
				iniS <- list()
				for (i in 1:length(chrOrder)){
					iniS[[i]] <- object@iniCall[chrDf$start[i]:chrDf$end[i]]
				}
				
				resFL <- lapply(iniS,cn.mops::segment,segMedianT=robust, minSeg=
								minWidth,...)
				
				for (i in 1:length(chrOrder)){
					resFL[[i]]$start <- resFL[[i]]$start+chrDf$start[i]-1
					resFL[[i]]$end <- resFL[[i]]$end+chrDf$start[i]-1
					resFL[[i]]$chr <- chrOrder[i]
				}
				
				segDf <- do.call(rbind,resFL)
				idx <- which(segDf$mean>robust)
				
				if (length(idx)==0){
					object@cnvr <- GRanges()
					message("No CNV regions detected with this settings.")
					return(object)	
				} else {
					segDf <- segDf[idx, ]
					cnvr <- GRanges(seqnames=segDf$chr,
							IRanges(start(object@gr)[segDf$start],
									end(object@gr)[segDf$end]))
					CN <- matrix(NA,nrow=length(idx),
							ncol(object@normalizedData))
					colnames(CN) <- colnames(object@normalizedData)
					GenomicRanges::values(cnvr) <- CN
					
					object@cnvr <- cnvr
					return(object)	
				}
			} else {
				if (length(object@cnvs)==0){
					object@cnvr <- GRanges()
					message("No individual CNVs in the result object.")
					return(object)	
				} else {
					cnvr <- GenomicRanges::reduce(object@cnvs)
					CN <- matrix(NA,nrow=length(cnvr),
							ncol(object@normalizedData))
					GenomicRanges::values(cnvr) <- CN
					object@cnvr <- cnvr
					return(object)
				}
			}
		})


