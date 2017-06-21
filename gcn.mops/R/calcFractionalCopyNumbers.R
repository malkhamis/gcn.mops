#' @title Calculation of fractional copy numbers for the CNVs and CNV regions.
#' 
#' @description This generic function calculates the fractional copy numbers of
#'  a CNV detection method stored in an instance of 
#' \code{\link{CNVDetectionResult-class}}. Must be a result of 
#' "referencecn.mops". 
#' 
#' @param object An instance of "CNVDetectionResult"
#' @param segStat Which statistic per segment should be used. Can be either
#' "mean" or "median". (Default="mean"). 
#' @examples
#' data(cn.mops)
#' r <- referencecn.mops(X[,1:2],apply(X,1,median))
#' calcFractionalCopyNumbers(r)
#' @return \code{calcFractionalCopyNumbers} returns an 
#' instance of "CNVDetectionResult".
#' @author Guenter Klambauer \email{klambauer@@bioinf.jku.at}
#' @export
#' @importFrom IRanges findOverlaps
#' @importFrom IRanges as.list
#' @importFrom IRanges as.matrix
#' @importFrom GenomicRanges values


setMethod("calcFractionalCopyNumbers", signature="CNVDetectionResult",
		definition = function(object,segStat="mean"){
			message(paste("This is an experimental function that has not",
							"been tested thoroughly! Please report bugs to",
							"the maintainer!"))
			
			
			priorImpact <- object@params$priorImpact
			cyc <- object@params$cyc
			classes <- object@params$classes
			I <- object@params$folds
			minReadCount <- object@params$minReadCount
			X <- object@normalizedData
			cnvr <- object@cnvr
			segmentation <- segmentation(object)
			cnvs <- cnvs(object)
			gr <- object@gr
			cov <- object@params$cov
			uT <- object@params$upperThreshold
			lT <- object@params$lowerThreshold
			mainClass <- object@params$mainClass
			method <- object@params$method
			mainCN <- 2*I[which(classes==mainClass)]
			
			usedMethod <- switch(method,
					"cn.mops"=.cn.mopsC,
					"haplocn.mops"=haplocn.mopsC,
					"referencecn.mops"=.referencecn.mops)
			segStatI <- switch(segStat,"mean"=mean,"median"=median)
			
			if (length(cnvr)==0 | length(cnvs)==0)
				stop(paste("No CNV regions in result object. Rerun cn.mops",
								"with different parameters!"))
			
			# for CNV regions
			M <- IRanges::as.list(IRanges::findOverlaps(cnvr,object@gr))
			XX <- lapply(M,function(i){ 
						if (length(i)>=3) ii <- i[-c(1,length(i))]
						else ii <- i
						apply(X[ii, ,drop=FALSE],2,segStatI) })
			if (method=="referencecn.mops"){
				lambda <- object@params$L[,mainClass]
				ll <- lapply(M,function(i){ 
							if (length(i)>=3) ii <- i[-c(1,length(i))]
							else ii <- i
							segStatI(lambda[ii]) 
						})
#				CN <-t(sapply(1:length(XX),function(j) {
#									round(mainCN*XX[[j]]/ll[[j]],1)								
#								}))
	alpha.prior <- rep(1,length(I))
	alpha.prior[which(classes=="CN2")] <- 1+priorImpact
	alpha.prior <- alpha.prior/sum(alpha.prior)
	
				CN <-t(sapply(1:length(XX),function(j) {
									paste("CN",format(
									round(mainCN*2^(usedMethod(x=XX[[j]],
													lambda=ll[[j]],
													I=I,
													classes=classes,
													cov=cov,
													minReadCount=minReadCount, 
													alpha.prior=alpha.prior,
													)$sini),1),nsmall=1),sep="")
								}))
				
			} else {
				stop("Only possible for referencecn.mops!")
			}
			
			CN <- matrix(CN,ncol=ncol(X))
			colnames(CN) <- colnames(X)
			resObject <- object
			GenomicRanges::values(cnvr) <- CN
			resObject@cnvr <- cnvr
			
			## now for individual CNVs and segmentation
			#mapping from CNVs to segmentation
			iCN <- rep(paste("CN",format(mainCN,nsmall=1),sep=""),length(segmentation))
			csM <- IRanges::as.matrix(IRanges::findOverlaps(segmentation,cnvs,type="within"))
			tmpIdx <- which(IRanges::values(segmentation)$sampleName[csM[,1]]==values(cnvs)$sampleName[csM[,2]])
			csM <- csM[tmpIdx, ,drop=FALSE]
			
			M2 <- IRanges::as.data.frame(segmentation)
			
			if (method=="referencecn.mops"){
				if (segStat=="mean"){
					CN2 <- paste("CN",format(round(2^(M2$mean)*mainCN,1),nsmall=1),
							sep="")
					
				} else if (segStat=="median"){
					CN2 <- paste("CN",format(round(2^(M2$median)*mainCN,1),nsmall=1),
							sep="")
					
				} else {
					stop("\"segStat\" must be \"mean\" or \"median\"!")
				}
			} else {
				stop("Only possible for referencecn.mops!")
			}
			
			
			GenomicRanges::values(segmentation)$CN <- CN2
			
			GenomicRanges::values(cnvs)$CN <- CN2[csM[,1]]
			resObject@cnvs <- cnvs
			resObject@segmentation <- segmentation
			
			return(resObject)							
		})
