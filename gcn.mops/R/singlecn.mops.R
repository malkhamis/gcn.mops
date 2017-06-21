# Copyright (C) 2013 Klambauer Guenter 
# <klambauer@bioinf.jku.at>

#cn.mops for given references
.singlecn.mops <- function(x,lambda,I = c(0.025,0.5,1,1.5,2,2.5,3,3.5,4), 
		classes=c("CN0","CN1","CN2","CN3","CN4","CN5","CN6","CN7","CN8"), cov,
		minReadCount=1,returnPosterior=FALSE,
		priorImpact=1) {
	N <- length(x)
	
	if (all(x<=minReadCount) & lambda <= minReadCount) {
		n <- length(I)
		idxCN2 <- which(classes=="CN2")
		alpha.est <- rep(0,n)
		alpha.est[idxCN2] <- 1
		expCN <- rep(classes[idxCN2],N)
		ini <- 0
		ExpLogFoldChange <- rep(0,N)
		if (returnPosterior){
			post.ik <- matrix(0,nrow=n,ncol=N)
			post.ik[idxCN2, ] <- 1
			l <-  list ("lambda"=lambda*I, "alpha"=alpha.est, "expectedCN"=expCN,
					"sini"=ExpLogFoldChange,"ini"=ini,"post"=post.ik)
			return(l) 
		} else {
			l <-  list ("lambda"=lambda*I, "alpha"=alpha.est, "expectedCN"=expCN,
					"sini"=ExpLogFoldChange,"ini"=ini,"post"=NA)
			return(l) 
		}
		
	} else {
		
		#browser()
		if (N>1){
			alpha.ik <- (sapply(I,function(ii) 
									exp(x*log(lambda*cov*ii)-lgamma(x+1)-lambda*cov*ii)
						))
			alpha.est <- rep(1,length(I))
			alpha.est[which(classes=="CN2")] <- 1+priorImpact
			alpha.est <- alpha.est/sum(alpha.est)
			alpha.ik <- pmax(alpha.ik,1e-15)
			alpha.ik <- alpha.ik*alpha.est
			
			alpha.ik <- t(alpha.ik/rowSums(alpha.ik))
			alpha.est <- rowMeans(alpha.ik)
			expCN <- classes[apply(alpha.ik,2,function(x) which(x==max(x))[1] )]
			ini <- mean(abs(log2(I)) %*% alpha.ik)
			ExpLogFoldChange <-  log2(I) %*%  alpha.ik
		} else {
			alpha.ik <- (sapply(I,function(ii) 
									exp(x*log(lambda*cov*ii)-lgamma(x+1)-lambda*cov*ii)
						))
			alpha.est <- rep(1,length(I))
			alpha.est[which(classes=="CN2")] <- 1+priorImpact
			alpha.est <- alpha.est/sum(alpha.est)
			alpha.ik <- pmax(alpha.ik,1e-15)
			alpha.ik <- alpha.ik*alpha.est
			alpha.ik <- alpha.ik/sum(alpha.ik)
			alpha.est <- (alpha.ik)
			expCN <- classes[ which.max(alpha.ik)][1]
			ini <- mean(abs(log2(I)) %*% alpha.ik)
			ExpLogFoldChange <-  log2(I) %*%  alpha.ik
		}	
		
		if (returnPosterior){
			l <-  list ("lambda"=lambda*I, "alpha"=alpha.est, "expectedCN"=expCN, 
					"sini"=ExpLogFoldChange, "ini"=ini, "post"=alpha.ik)
		} else {
			l <-  list ("lambda"=lambda*I, "alpha"=alpha.est, "expectedCN"=expCN, 
					"sini"=ExpLogFoldChange, "ini"=ini, "post"=NA)
		}
		return(l)
	}
}




#' @title Copy number detection in NGS data with in a setting in which only 
#' one sample is available
#' @name singlecn.mops
#' 
#' @description This function performs the an alternative version of the
#' cn.mops algorithm adapted to a setting of a single sample.
#' 
#' @param x Either an instance of "GRanges" or a raw data matrix with one column
#' or a vector of read counts. An entry is the read count of the  sample 
#' in the genomic region.
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
#' be assumed to have copy number 2. Default = 1.
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
#' for copy number 3 or 4. Default = 0.5.
#' @param lowerThreshold Negative real value that sets the cut-off for copy
#' number losses. All CNV calling values below this value will be called as 
#' "loss". The value should be set close to the log2 of the expected foldchange
#' for copy number 1 or 0. Default = -0.9.
#' @param minWidth Positive integer that is exactly the parameter "min.width"
#' of the "segment" function of "DNAcopy". minWidth is the minimum number 
#' of segments a CNV should span. Default = 3.
#' @param segAlgorithm Which segmentation algorithm should be used. If set to
#' "DNAcopy" circular binary segmentation is performed. Any other value will
#' initiate the use of our fast segmentation algorithm. Default = "fast".
#' @param minReadCount If all samples are below this value the algorithm will
#' return the prior knowledge. This prevents that the algorithm from being 
#' applied to segments with very low coverage. Default=1. 
#' @param returnPosterior Flag that decides whether the posterior probabilities
#' should be returned. The posterior probabilities have a dimension of samples
#' times copy number states times genomic regions and therefore consume a lot
#' of memory. Default=FALSE.
#' @param ... Additional parameters will be passed to the "DNAcopy"
#' or the standard segmentation algorithm.
#' @examples 
#' data(cn.mops)
#' singlecn.mops(XRanges[,1])
#' @importFrom parallel makeCluster
#' @importFrom parallel clusterEvalQ
#' @importFrom parallel parApply
#' @importFrom parallel stopCluster 
#' @return An instance of "CNVDetectionResult".
#' @author Guenter Klambauer \email{klambauer@@bioinf.jku.at}
#' @export

singlecn.mops <- function(x,I = c(0.025,0.5,1,1.5,2,2.5,3,3.5,4),
		classes=c("CN0","CN1","CN2","CN3","CN4","CN5","CN6","CN7","CN8"),
		priorImpact = 1,cyc = 20,parallel=0,
		norm=1, normType="poisson",sizeFactor="mean",normQu=0.25, quSizeFactor=0.75,
		upperThreshold=0.5,lowerThreshold=-0.9,
		minWidth=3,segAlgorithm="fast",minReadCount=1,
		returnPosterior=FALSE,...){
	
	
	############ check input ##################################################
	if (is.vector(x)){
		inputType <- "DataMatrix"
		X <- matrix(x,ncol=1)
		chr <- rep("undef",nrow(X))
		irAllRegions <- IRanges(start=1:nrow(X),end=1:nrow(X))
		grAllRegions <- GRanges(chr,irAllRegions)
		start <- 1:nrow(X)
		end <- 1:nrow(X)
		rownames(X) <- paste(chr,start,end,sep="_")
		
	} else if(class(x)=="GRanges"){
	
		inputType <- "GRanges"
		x <- sortSeqlevels(x)
		X <- IRanges::as.matrix(IRanges::values(x))
		
		if (ncol(as.matrix(values(x)))!=1){
			stop("More than one sample detected. Use standard cn.mops!")
		}
		if (length(IRanges::unique(strand(x))) >1){
			stop(paste("Different strands found in GRanges object. Please make",
							"read counts independent of strand."))
		}
		chr <- as.character(seqnames(x))
		start <- start(x)
		end <- end(x)
		rownames(X) <- paste(chr,start,end,sep="_")
		
		irAllRegions <- IRanges(start,end)
		grAllRegions <- GRanges(chr,irAllRegions,seqinfo=seqinfo(x))
		grAllRegions <- sortSeqlevels(grAllRegions) 
		names(irAllRegions) <- NULL
		
	} else if (is.matrix(x)){
		if (nrow(x)!=1){
			stop("More than one sample detected. Use standard cn.mops!")
		}
		
		inputType <- "DataMatrix"
		X <- x
		chr <- rep("undef",nrow(X))
		irAllRegions <- IRanges(start=1:nrow(X),end=1:nrow(X))
		grAllRegions <- GRanges(chr,irAllRegions)
		start <- 1:nrow(X)
		end <- 1:nrow(X)
		rownames(X) <- paste(chr,start,end,sep="_")
		
	} else{
		stop("GRanges object or read count matrix needed as input.")
	}
	
	
	if (any(X<0) | any(!is.finite(X))){
		stop("All values must be greater or equal zero and finite.\n")
	}
	
	
	if (!is.numeric(I)) stop("\"I\" must be numeric.")
	if (!is.character(classes)) stop("\"classes\" must be character.")
	if (length(I)!=length(classes)){
		stop("I and classes must have same length!")
	}
	if (!("CN2" %in% classes)){stop("One element of classes must be CN2 .\n")}
	if (!(is.numeric(priorImpact) & length(priorImpact)==1)) 
		stop("\"priorImpact\" be must numeric and of length 1.")
	if (!(is.numeric(cyc) & length(cyc)==1)) 
		stop("\"cyc\" must be numeric and of length 1.")
	if (!(is.numeric(parallel) & length(parallel)==1)) 
		stop("\"parallel\" must be numeric and of length 1.")
	if (!(normType %in% c("mean","min","median","quant","mode","poisson"))){
		stop(paste("Set normalization to \"mean\"",
						"\"min\", \"median\", \"quant\" or \"mode\"."))
	}
	if (!(is.numeric(normQu) & length(normQu)==1)) 
		stop("\"normQu\" must be numeric and of length 1.")
	if (is.logical(norm))
		norm <- as.integer(norm)
	if (!(norm %in% c(0,1,2)))
		stop("\"norm\" must be 0,1 or 2.")
	if (!(is.numeric(lowerThreshold) & length(lowerThreshold)==1)) 
		stop("\"lowerThreshold\" must be numeric and of length 1.")
	if (!(is.numeric(upperThreshold) & length(upperThreshold)==1)) 
		stop("\"upperThreshold\" must be numeric and of length 1.")
	if (!(is.numeric(minWidth) & length(minWidth)==1)) 
		stop("\"minWidth\" must be numeric and of length 1.")
	if (!is.character(segAlgorithm)){
		stop("\"segAlgorithm\" must be \"fastseg\" or \"DNAcopy\"!")
	}
	if (!(is.numeric(minReadCount) & length(minReadCount)==1)) 
		stop("\"minReadCount\" must be numeric and of length 1.")
	if (!is.logical(returnPosterior))
		stop("\"returnPosterior\" must be logical.")	
	
	
	if (is.null(colnames(X))){
		colnames(X) <- paste("Sample",1:ncol(X),sep="_")
	}
	############################################################################
	
	version <- packageDescription("cn.mops")$Version
	params <- list("referencecn.mops",I,
			classes,
			priorImpact,cyc,
			normType,normQu,
			upperThreshold,lowerThreshold,
			minWidth,segAlgorithm,minReadCount,"CN2",version,paste(...))
	names(params) <- c("method","folds",
			"classes",
			"priorImpact","cyc",
			"normType","normQu",
			"upperThreshold","lowerThreshold",
			"minWidth","segAlgorithm","minReadCount","mainClass",
			"version","SegmentationParams")
	############################################################################
	m <- nrow(X)
	N <- ncol(X)
	
	n <- length(I)
	chrOrder <- unique(chr) #unique chromosome names in alphabetic order
	chrBpts <- cumsum(table(chr)[chrOrder])
# contains the idx where chromosomes start and end in X
	chrDf <- data.frame(start=c(1,chrBpts[-length(chrBpts)]+1),
			end=chrBpts)
	rownames(chrDf) <- chrOrder
	
	message(paste("Normalization is not  applicable for singlecn.mops! \n",
						"Normalization parameters are only dummies."))
	
	
	if (norm==0){
		X.norm <- X
		cov <- rep(1,N)
	} else if (norm==1) {
		message("Normalizing...")
		X.norm <- X
		
		cov <- rep(1,N)
	} else if (norm==2) {
		X.viz <- X		
		X.norm <- X
		cov <- rep(1,N)
		
	} else {
		stop("\"norm\" must be 0,1 or 2.")
	}
	params$cov <- cov
	
	message("Starting local modeling, please be patient...  ")
	
	if (parallel > 0){
		cl <- parallel::makeCluster(as.integer(parallel),type="SOCK")
		parallel::clusterEvalQ(cl,".referencecn.mops")
	} 
	
	res <- list()
	#browser()
	#lambda <- rep(median(X[which(X[,1]>0),1]),m)
	lambda <- unlist(apply(chrDf,1,function(ii){
				yy <- X[ii[1]:ii[2], ]
				mm <- median(yy[which(yy>0)])
				return(rep(mm,ii[2]-ii[1]+1))
			}))
	
	
	for (chrom in chrOrder){
		message(paste("Reference sequence: ",chrom))
		chrIdx <- chrDf[chrom,1]:chrDf[chrom,2]
		
		#cn.mops params
		#n,N, I set
		idxCN2 <- which(classes=="CN2")
		
		#function(x,lambda,I = c(0.025,0.5,1,1.5,2,2.5,3,3.5,4), 
		#		classes=c("CN0","CN1","CN2","CN3","CN4","CN5","CN6","CN7","CN8"), cov,
		#			minReadCount=1,returnPosterior=FALSE)
		
		if (parallel==0){
			resChr <- lapply(1:length(chrIdx),function(i)
						.singlecn.mops(X.norm[chrIdx[i], ,drop=FALSE],lambda=
										lambda[chrIdx[i]],I=I,classes=classes,cov=cov,
								minReadCount=minReadCount,
								returnPosterior=returnPosterior,
								priorImpact=priorImpact))
		} else {
			
			resChr <- parallel::parLapply(cl,1:length(chrIdx),function(i)
						.singlecn.mops(X.norm[chrIdx[i], ,drop=FALSE],lambda=
										lambda[chrIdx[i]],I=I,classes=classes,cov=cov,
								minReadCount=minReadCount,
								returnPosterior=returnPosterior,
								priorImpact=priorImpact))
		}
		
		res <- c(res, resChr)
	}
	if (parallel > 0){
		parallel::stopCluster(cl)
	} 
	
	## Postprocess result
	L <- t(sapply(res,.subset2,1))
	rownames(L) <- rownames(X)
	colnames(L) <- classes
	params$L <- L
	A <- t(sapply(res,.subset2,2))
	rownames(A) <- rownames(X)
	colnames(A) <- classes
	params$A <- A
	if (N==1){
		CN <- matrix(sapply(res,.subset2,3),ncol=1)
		rownames(CN) <- rownames(X)
		colnames(CN) <- colnames(X)
		sINI <- matrix(sapply(res,.subset2,4),ncol=1)
		rownames(sINI) <- rownames(X)
		colnames(sINI) <- colnames(X)
	} else {
		CN <- t(sapply(res,.subset2,3))
		rownames(CN) <- rownames(X)
		colnames(CN) <- colnames(X)
		sINI <- t(sapply(res,.subset2,4))
		rownames(sINI) <- rownames(X)
		colnames(sINI) <- colnames(X)
	}
	INI <- (sapply(res,.subset2,5))
	names(INI) <- rownames(X)
	
	if (returnPosterior){
		tt <- try(post <- array(dim=c(m,n,N)))
		if (inherits(tt,"try-error")){
			message("Posterior too large for array extent.")
			post <- array(NA,dim=c(1,1,1))
		} else {
			post.tmp <- t(lapply(res,.subset2,6))
			for (i in 1:m){
				post[i, ,] <- post.tmp[[i]]
			}
			dimnames(post) <- list(NULL,classes,colnames(X))
			rm("post.tmp")
		}
	} else {
		post <- array(NA,dim=c(1,1,1))
	}
	rm("res")
	
	
	if (m>5){
		message("Starting segmentation algorithm...")
		
		if (segAlgorithm=="DNAcopy"){
			message("Using \"DNAcopy\" for segmentation.")
			requireNamespace("DNAcopy")
			if (!exists("eta")){eta <- 0.05}
			if (!exists("nperm")){nperm <- 10000}
			if (!exists("alpha")){alpha <- 0.01}
			if (minWidth > 5){
				message("For DNAcopy the maximum \"minWidth\" is 5.")
				message("Resetting \"minWidth\" to 5.")
				minWidth <- 5
			}
			DNAcopyBdry <- DNAcopy::getbdry(eta=eta,nperm=nperm,tol=alpha,
					max.ones=floor(nperm*alpha)+1)
			
			
			if (parallel==0){
				resSegm <- apply(sINI,2,.segmentation,
						chr=chr,minWidth=minWidth,DNAcopyBdry=DNAcopyBdry,...)
			} else {
				cl <- parallel::makeCluster(as.integer(parallel),type="SOCK")
				parallel::clusterEvalQ(cl,".segmentation")
				resSegm <- parallel::parApply(cl,sINI,2,.segmentation,
						chr=chr,minWidth=minWidth,DNAcopyBdry=DNAcopyBdry,...)
				parallel::stopCluster(cl)
			}
			
			
			
			#resSegm <- lapply(resSegm,function(x) x <- x[order(x$chr,x$start), ])
			segDf <- cbind(do.call(rbind,resSegm),
					rep(colnames(X),sapply(resSegm,nrow)))
			rm("resSegm")
			
			segDf$CN <- NA
			
			colnames(segDf) <-
					c("chr","start","end","mean","median","sample","CN")
			segDf <- segDf[ ,c("chr","start","end","sample","median","mean","CN")]
			segDf <- segDf[order(as.character(segDf$chr),segDf$sample,segDf$start), ]
			
			
			callsS <- matrix(NA,nrow=m,ncol=N)
			for (chrom in chrOrder){
				chrIdx <- chrDf[chrom,1]:chrDf[chrom,2]
				segDfTmp <- subset(segDf,chr==chrom)
				callsS[chrIdx, ] <- 
						matrix(rep(segDfTmp$mean,segDfTmp$end-segDfTmp$start+1),
								ncol=N)
			}
			
			colnames(callsS) <- colnames(X)
			
			
			segDfSubset <- segDf[which(
							segDf$mean >= upperThreshold |
									segDf$mean <= lowerThreshold), ]
			segDfSubset <- segDfSubset[which(
							(segDfSubset$end-segDfSubset$start+1) >= minWidth), ]
			
			
			
		} else {
			message("Using \"fastseg\" for segmentation.")
			resSegmList <- list()
			segDf <- data.frame(stringsAsFactors=FALSE)
			
			callsS <- matrix(NA,nrow=m,ncol=N)
			colnames(callsS) <- colnames(X)
			for (chrom in chrOrder){
				chrIdx <- chrDf[chrom,1]:chrDf[chrom,2]
				
				if (parallel==0){
					resSegmList[[chrom]] <- apply(sINI[chrIdx, ,drop=FALSE],2,
							segment,
							minSeg=minWidth,...)
				} else {
					cl <- parallel::makeCluster(as.integer(parallel),type="SOCK")
					parallel::clusterEvalQ(cl,"segment")
					resSegmList[[chrom]] <- parallel::parApply(cl,sINI[chrIdx, ,drop=FALSE],2,
							segment,minSeg=minWidth,...)
					parallel::stopCluster(cl)
				}
				
				segDfTmp <- cbind(do.call(rbind,resSegmList[[chrom]]),
						"sample"=rep(colnames(X),
								sapply(resSegmList[[chrom]],nrow)))
				segDfTmp$chr <- chrom
				
				
				callsS[chrIdx, ] <- 
						matrix(rep(segDfTmp$mean,segDfTmp$end-segDfTmp$start+1),
								ncol=N)
				
				segDf <- rbind(segDf,segDfTmp)
			}
			segDf <- data.frame(segDf,"CN"=NA,stringsAsFactors=FALSE)		
			colnames(segDf) <- c("start","end","mean","median","sample",
					"chr","CN")
			
			segDfSubset <- segDf[which(segDf$mean >= upperThreshold
									| segDf$mean <= lowerThreshold), ]	
			segDfSubset <- 
					segDfSubset[which((segDfSubset$end-segDfSubset$start+1)
											>= minWidth), ]
		}
		
		
		if (nrow(segDfSubset)>0){
			
			
			# Assembly of result object
			r <- new("CNVDetectionResult")
			
			sampleNames <- segDfSubset$sample
			
			if (inputType=="GRanges"){
				ir <- IRanges()
				for (chrom in chrOrder){
					chrIdx <- chrDf[chrom,1]:chrDf[chrom,2]
					inputChr <- x[chrIdx]
					segDfSubsetChr <- subset(segDfSubset,chr==chrom)
					if (nrow(segDfSubsetChr) >0){
						ir <- c(ir,IRanges(start(inputChr)[
												segDfSubsetChr$start],
										end(inputChr)[segDfSubsetChr$end]))
					}
				}
			} else if (inputType=="DataMatrix"){
				ir <- IRanges(start=segDfSubset$start,end=segDfSubset$end)
			}
			
			
			rd <- GRanges(seqnames=segDfSubset$chr,ir,"sampleName"=sampleNames,
					"median"=segDfSubset$median,"mean"=segDfSubset$mean,
					"CN"=segDfSubset$CN,seqinfo=seqinfo(grAllRegions))
			rd <- sortSeqlevels(rd) 
			
			
#			cnvr <- GRanges(seqnames=seqnames(cnvrR),irCNVR)
#			GenomicRanges::values(cnvr) <- cnvrCN
			
			if (norm==2){
				r@normalizedData    <- X.viz		
			} else {
				r@normalizedData    <- X.norm			
			}
			r@localAssessments  <- sINI
			r@individualCall   	<- callsS
			r@iniCall        	<- INI
			r@cnvs				<- rd
			r@cnvr				<- rd
			
			
			
			if (inputType=="GRanges"){
				irS <- IRanges()
				for (chrom in chrOrder){
					chrIdx <- chrDf[chrom,1]:chrDf[chrom,2]
					inputChr <- x[chrIdx]
					segDfChr <- subset(segDf,chr==chrom)
					if (nrow(segDfChr) >0 ){
						irS <- c(irS,IRanges(start(inputChr)[segDfChr$start],
										end(inputChr)[segDfChr$end]))
					}
				}
				r@segmentation 			<- GRanges(seqnames=segDf$chr,
						irS,
						"sampleName"=segDf$sample,"median"=segDf$median,
						"mean"=segDf$mean,"CN"=segDf$CN,seqinfo=seqinfo(grAllRegions))
				r@segmentation <- sortSeqlevels(r@segmentation) 
				
			} else if (inputType=="DataMatrix"){
				r@segmentation <- GRanges(seqnames=segDf$chr,
						IRanges(segDf$start,segDf$end),"sampleName"=segDf$sample,
						"median"=segDf$median,
						"mean"=segDf$mean,"CN"=segDf$CN,seqinfo=seqinfo(grAllRegions))
				r@segmentation <- sortSeqlevels(r@segmentation) 
				
			}
			
			r@gr <- grAllRegions
			r@posteriorProbs 	<- post
			r@params			<- params
			r@integerCopyNumber	<- CN
			r@sampleNames		<- colnames(X)
			
			return(r)	
		} else {
			message(paste("No CNVs detected. Try changing \"normalization\",", 
							"\"priorImpact\" or \"thresholds\"."))
			# Assembly of result object
			r <- new("CNVDetectionResult")
			if (inputType=="GRanges"){
				irS <- IRanges()
				for (chrom in chrOrder){
					chrIdx <- chrDf[chrom,1]:chrDf[chrom,2]
					inputChr <- x[chrIdx]
					segDfChr <- subset(segDf,chr==chrom)
					if (nrow(segDfChr) >0 ){
						irS <- c(irS,IRanges(start(inputChr)[segDfChr$start],
										end(inputChr)[segDfChr$end]))
					}
				}
				r@segmentation 			<- GRanges(seqnames=segDf$chr,
						irS,
						"sampleName"=segDf$sample,"median"=segDf$median,
						"mean"=segDf$mean,"CN"=segDf$CN,seqinfo=seqinfo(grAllRegions))
				r@segmentation <- sortSeqlevels(r@segmentation) 
				
			} else if (inputType=="DataMatrix"){
				r@segmentation 			<- GRanges(seqnames=segDf$chr,
						IRanges(segDf$start,segDf$end),
						"sampleName"=segDf$sample,"median"=segDf$median,
						"mean"=segDf$mean,"CN"=segDf$CN,seqinfo=seqinfo(grAllRegions))
				r@segmentation <- sortSeqlevels(r@segmentation) 
				
			}
			
			r@gr <- grAllRegions
			if (norm==2){
				r@normalizedData    <- X.viz		
			} else {
				r@normalizedData    <- X.norm			
			}
			r@localAssessments  <- sINI
			r@individualCall   	<- callsS
			r@iniCall        	<- INI
			r@posteriorProbs 	<- post
			r@params			<- params
			r@integerCopyNumber	<- CN
			r@sampleNames		<- colnames(X)
			return(r)	
		}
		
		
	} else {
		message(paste("Less than five genomic segments considered, therefore no",
						" segmentation."))	
		# Assembly of result object
		r <- new("CNVDetectionResult")	#
		r@gr <- grAllRegions
		if (norm==2){
			r@normalizedData    <- X.viz		
		} else {
			r@normalizedData    <- X.norm			
		}
		r@localAssessments  <- sINI
		r@individualCall   	<- sINI
		r@params			<- params
		r@integerCopyNumber	<- CN
		r@sampleNames		<- colnames(X)
		return(r)	
		
	}
	
	
	message("GC correction is not implemented yet and will improve results.")
	
}

