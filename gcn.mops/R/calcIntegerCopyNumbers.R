#' @title Calculation of integer copy numbers for the CNVs and CNV regions.
#' @description This generic function calculates the integer copy numbers of
#'  a CNV detection method stored in an instance of 
#' \code{\link{CNVDetectionResult-class}}. 
#' 
#' @param object An instance of "CNVDetectionResult" 
#' @examples
#' data(cn.mops)
#' r <- cn.mops(X[1:100,1:5])
#' calcIntegerCopyNumbers(r)
#' @return \code{calcIntegerCopyNumbers} returns an 
#' instance of "CNVDetectionResult".
#' @author Guenter Klambauer \email{klambauer@@bioinf.jku.at}
#' @export
#' @importFrom IRanges findOverlaps
#' @importFrom IRanges as.list
#' @importFrom IRanges as.matrix
#' @importFrom GenomicRanges values


setMethod("calcIntegerCopyNumbers", signature="CNVDetectionResult",
		definition = function(object){
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
			mainClassIdx <- which(classes==mainClass)
			method <- object@params$method
			if ("useMedian" %in% names(object@params)){
				useMedian <- object@params$useMedian				
			} else {
				useMedian <- FALSE
			}
			
			usedMethod <- switch(method,
					"cn.mops"=.cn.mopsC,
					"haplocn.mops"=haplocn.mopsC,
					"referencecn.mops"=.referencecn.mops)
			
			
			if (length(cnvr)==0 | length(cnvs)==0)
				stop(paste("No CNV regions in result object. Rerun cn.mops",
								"with different parameters!"))
			
			### now for individual CNVs and segmentation			
			iCN <- rep(mainClass,length(segmentation))
			#mapping from CNVs to segmentation
			csM <- IRanges::as.matrix(IRanges::findOverlaps(segmentation,cnvs,type="within"))
			
			
			tmpIdx <- which(values(segmentation)$sampleName[csM[,1]]==values(cnvs)$sampleName[csM[,2]])
			csM <- csM[tmpIdx, ,drop=FALSE]
			
			# check again for equality
			csM <- csM[(segmentation[csM[,1]]==cnvs[csM[,2]]), ,drop=FALSE]
			
			idx <- csM[,1]
			
			M2 <- IRanges::as.list(IRanges::findOverlaps(segmentation[idx],object@gr))
			XX2 <- lapply(M2,function(i){ 
						if (length(i)>=3) ii <- i[-c(1,length(i))]
						else ii <- i
						apply(X[ii, ,drop=FALSE],2,mean) })
			
			if (method=="referencecn.mops"){
				lambda <- object@params$L[,mainClass]
				ll <- lapply(M2,function(i){ 
							if (length(i)>=3) ii <- i[-c(1,length(i))]
							else ii <- i
							mean(lambda[ii]) 
						})
				alpha.prior <- rep(1,length(I))
				alpha.prior[which(classes=="CN2")] <- 1+priorImpact
				alpha.prior <- alpha.prior/sum(alpha.prior)
				
				CN2 <-t(sapply(1:length(XX2),function(j) {
									usedMethod(x=XX2[[j]],
											lambda=ll[[j]],
											I=I,
											classes=classes,
											cov=cov,
											minReadCount=minReadCount,
											alpha.prior=alpha.prior)$expectedCN
								}))
			} else {
				CN2 <-t(sapply(XX2,function(x) usedMethod(x,I=I,
											classes=classes,
											cov=cov,priorImpact=priorImpact,
											cyc=cyc,
											minReadCount=minReadCount)$expectedCN))
			}
			
			CN2 <- matrix(CN2,ncol=ncol(X))
			colnames(CN2) <- colnames(X)
			extractedCN <- CN2[cbind(1:length(idx),
							match(as.character(values(segmentation[idx])$sampleName),
									colnames(X)))]
			
			iCN[idx] <- extractedCN 
			GenomicRanges::values(segmentation)$CN <- iCN
			
			
			if (useMedian){		
				GenomicRanges::values(segmentation)$CN[which(values(segmentation)$CN==mainClass & values(segmentation)$median > uT)] <-
						classes[mainClassIdx+1]
				GenomicRanges::values(segmentation)$CN[which(values(segmentation)$CN==mainClass & values(segmentation)$median < lT)] <- 
						classes[mainClassIdx-1]
				
				GenomicRanges::values(cnvs)$CN <- extractedCN[csM[,2]]
				GenomicRanges::values(cnvs)$CN[which(values(cnvs)$CN==mainClass & values(cnvs)$median > uT)] <-
						classes[mainClassIdx+1]
				GenomicRanges::values(cnvs)$CN[which(values(cnvs)$CN==mainClass & values(cnvs)$median < lT)] <- 
						classes[mainClassIdx-1]
				
			}  else {
				GenomicRanges::values(segmentation)$CN[which(values(segmentation)$CN==mainClass & values(segmentation)$mean > uT)] <-
						classes[mainClassIdx+1]
				GenomicRanges::values(segmentation)$CN[which(values(segmentation)$CN==mainClass & values(segmentation)$mean < lT)] <- 
						classes[mainClassIdx-1]
				
				GenomicRanges::values(cnvs)$CN <- extractedCN[csM[,2]]
				GenomicRanges::values(cnvs)$CN[which(values(cnvs)$CN==mainClass & values(cnvs)$mean > uT)] <-
						classes[mainClassIdx+1]
				GenomicRanges::values(cnvs)$CN[which(values(cnvs)$CN==mainClass & values(cnvs)$mean < lT)] <- 
						classes[mainClassIdx-1]
				
			}
			
			### for CNV regions
			
			
			M <- IRanges::as.list(IRanges::findOverlaps(cnvr,cnvs))
			CN <-t(sapply(1:length(M),function(i){
								xxx <- rep(mainClass,ncol(X))
								names(xxx) <- colnames(X)
								zz <- cnvs[M[[i]]]$CN
								zz2 <- as.character(cnvs[M[[i]]]$sampleName)
								zIdx <- match(unique(zz2),zz2)
								xxx[zz2[zIdx]] <- zz[zIdx] 
								return(xxx)
							}))
			
			
			if (ncol(X)==1)
				CN <- matrix(CN,ncol=1)
			colnames(CN) <- colnames(X)
			
#			XX <- lapply(M,function(i){ 
#						if (length(i)>=3) ii <- i[-c(1,length(i))]
#						else ii <- i
#						apply(X[ii, ,drop=FALSE],2,mean) })
#			if (method=="referencecn.mops"){
#				lambda <- object@params$L[,mainClass]
#				ll <- lapply(M,function(i){ 
#							if (length(i)>=3) ii <- i[-c(1,length(i))]
#							else ii <- i
#							mean(lambda[ii]) 
#						})
#				CN <-t(sapply(1:length(XX),function(j) {
#								usedMethod(x=XX[[j]],
#								lambda=ll[[j]],
#								I=I,
#								classes=classes,
#								cov=cov,
#								minReadCount=minReadCount)$expectedCN
#								}))
#				
#			} else {
#				CN <-t(sapply(XX,function(x) usedMethod(x,I=I,
#											classes=classes,
#											cov=cov,priorImpact=priorImpact,
#											cyc=cyc,
#											minReadCount=minReadCount)$expectedCN))
#			}
#			
#			CN <- matrix(CN,ncol=ncol(X))
#			colnames(CN) <- colnames(X)
			
			
			
			
			resObject <- object
			GenomicRanges::values(cnvr) <- CN
			resObject@cnvr <- cnvr
			resObject@cnvs <- cnvs
			resObject@segmentation <- segmentation
			return(resObject)							
		})

