# Copyright (C) 2011 Klambauer Guenter 
# <klambauer@bioinf.jku.at>

haplocn.mopsC <- function(x,I = c(0.025,1,2,3,4,5,6,7,8), 
		classes=c("CN0","CN1","CN2","CN3","CN4","CN5","CN6","CN7","CN8"), cov,
		priorImpact = 1,cyc = 20, minReadCount=1) {
	
	
	N <- length(x)
	n <- length(I)
	
	if (missing(cov)){cov <- rep(1,N)}
	
	idxCN1 <- which(classes=="CN1")
	
	alpha.init <- rep(0.05,n)
	alpha.init[idxCN1] <- 0.6
	alpha.init <- alpha.init/ sum(alpha.init)
	alpha.est <- alpha.init
	alpha.prior <- rep(0,n)
	alpha.prior[idxCN1] <- 1
	alpha.prior <- alpha.prior*priorImpact
	
	if (all(x<=minReadCount)) {
		lambda.est <- rep(0,n)
		alpha.est <- rep(0,n)
		alpha.est[idxCN1] <- 1
		expCN <- rep(classes[idxCN1],N)
		ini <- 0
		ExpLogFoldChange <- rep(0,N)
		post.ik <- matrix(0,nrow=n,ncol=N)
		post.ik[idxCN1, ] <- 1
		
		
		params <- list(n,classes,I,priorImpact,cyc)
		names(params) <- c("nclasses","classes","I","priorImpact","cyc")
		l <-  list ("lambda"=lambda.est, "alpha"=alpha.est, "expectedCN"=expCN,
				"sini"=ExpLogFoldChange,"ini"=ini,"post"=post.ik, 
				"params"=params)
		return(l)
		
	} else {
		lambda.est <- median(x*1/cov,na.rm=TRUE)
		if (lambda.est < 1e-10){lambda.est <- max(mean(x*1/cov,na.rm=TRUE),1.0)}
		lambda.init <- I*lambda.est
		ret=.Call("cnmops", as.numeric(x), I,cov, as.integer(cyc), alpha.init, lambda.init,
				alpha.prior)
		alpha.ik=ret$alpha.ik
		alpha.i=ret$alpha.i
		alpha.est=ret$alpha.est
		lambda.est=ret$lambda.est
		
		#Posterior
		post.ik <- alpha.ik
		if (is.null(names(x))){
			colnames(post.ik) <- paste("x_",1:N,sep="")
		} else {
			colnames(post.ik) <- names(x)
		}
		rownames(post.ik) <- classes
		expCN <- classes[apply(post.ik,2,function(x) which(x==max(x))[1] )]
		
		ini <- mean(abs(log2(I)) %*% alpha.ik)
		ExpLogFoldChange <-  log2(I) %*%  post.ik
		params <- list(n,classes,I,priorImpact,cyc)
		names(params) <- c("nclasses","classes","I","priorImpact","cyc")
		l <-  list ("lambda"=lambda.est, "alpha"=alpha.est, "expectedCN"=expCN, 
				"sini"=ExpLogFoldChange, "ini"=ini, "post"=post.ik,
				"params"=params,"version"=version)
		return(l)
	}
}

haplocn.mopsCE <- function(x, I, classes, cov, cyc, N, n,
		idxCN1, alphaInit, alphaPrior, minReadCount) {
	if (all(x<=minReadCount)) {
		lambda.est <- rep(0,n)
		alpha.est <- rep(0,n)
		alpha.est[idxCN1] <- 1
		expCN <- rep(classes[idxCN1],N)
		ini <- 0
		ExpLogFoldChange <- rep(0,N)
		post.ik <- matrix(0,nrow=n,ncol=N)
		post.ik[idxCN1, ] <- 1
		l <-  list ("lambda"=lambda.est, "alpha"=alpha.est, "expectedCN"=expCN,
				"sini"=ExpLogFoldChange,"ini"=ini,"post"=post.ik)
		return(l)
	} else {
		lambda.est <- median(x*1/cov,na.rm=TRUE)
		if (lambda.est < 1e-10){lambda.est <- max(mean(x*1/cov,na.rm=TRUE),1.0)}
		lambda.init <- I*lambda.est
		ret <- .Call("cnmops", as.numeric(x), I,cov, 
				as.integer(cyc), alphaInit, lambda.init,
				alphaPrior)
		expCN <- classes[apply(ret$alpha.ik,2,function(x) which(x==max(x))[1] )]
		ini <- mean(abs(log2(I)) %*% ret$alpha.ik)
		ExpLogFoldChange <-  log2(I) %*%  ret$alpha.ik
		l <-  list ("lambda"=ret$lambda.est, "alpha"=ret$alpha.est, "expectedCN"=expCN, 
				"sini"=ExpLogFoldChange, "ini"=ini, "post"=ret$alpha.ik)
		return(l)
	}
}


#' @title Copy number detection in NGS data of haploid samples.
#' 
#' @description Performs the cn.mops algorithm for copy number detection in
#' NGS data adjusted to haploid genomes. It is assumed that the normal state
#' is copy number 1. This is an experimental method at the moment.
#' 
#' @param input Either an instance of "GRanges" or a raw data matrix, where
#' columns are interpreted as samples and rows as genomic regions. An entry is
#' the read count of a sample in the genomic region.
#' @param I Vector positive real values that contain the expected fold change
#' of the copy number classes.  Length of this vector must be equal to the 
#' length of the "classes" parameter vector. For copy number polymorphisms 
#' in haploid organisms we suggest to use the default 
#' I = c(0.025,1,2,3,4,5,6,7,8).
#' @param classes Vector of characters of the same length as the parameter
#' vector "I". One vector element must be named "CN1". The names reflect the 
#' labels of the copy number classes. 
#' Default = c("CN0","CN1","CN2","CN3","CN4","CN5","CN6","CN7","CN8").
#' @param priorImpact Positive real value that reflects how strong the prior
#' assumption affects the result. The higher the value the more samples will
#' be assumed to have copy number 1. Default = 1.
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
#' of segments a CNV should span. Default = 4.
#' @param segAlgorithm Which segmentation algorithm should be used. If set to
#' "DNAcopy" circular binary segmentation is performed. Any other value will
#' initiate the use of our fast segmentation algorithm. Default = "fast".
#' @param minReadCount If all samples are below this value the algorithm will
#' return the prior knowledge. This prevents that the algorithm from being 
#' applied to segments with very low coverage.
#' @param ... Additional parameters will be passed to the "DNAcopy"
#' or the standard segmentation algorithm.
#' @param returnPosterior Flag that decides whether the posterior probabilities
#' should be returned. The posterior probabilities have a dimension of samples
#' times copy number states times genomic regions and therefore consume a lot
#' of memory. Default=FALSE.
#' @examples 
#' data(cn.mops)
#' haplocn.mops(XRanges[1:200, ])
#' haplocn.mops(XRanges[1:200, ],parallel=2)
#'
#' @importFrom parallel makeCluster
#' @importFrom parallel clusterEvalQ
#' @importFrom parallel parApply
#' @importFrom parallel stopCluster 
#' @useDynLib cn.mops
#' @return An instance of "CNVDetectionResult".
#' @author Guenter Klambauer \email{klambauer@@bioinf.jku.at}
#' @export


haplocn.mops <- function(input,I = c(0.025,1,2,3,4,5,6,7,8),
		classes=c("CN0","CN1","CN2","CN3","CN4","CN5","CN6","CN7","CN8"),
		priorImpact = 1,cyc = 20,parallel=0,
		norm=1, normType="poisson",sizeFactor="mean",normQu=0.25, quSizeFactor=0.75,
		upperThreshold=0.6,lowerThreshold=-0.9,
		minWidth=3,segAlgorithm="fast",minReadCount=1,
		returnPosterior=FALSE,...){
	
	version <- packageDescription("cn.mops")$Version
	
	
	############ check input ##################################################
	if(class(input)=="GRanges"){
		inputType <- "GRanges"
		input <- sortSeqlevels(input)
		X <- IRanges::as.matrix(IRanges::values(input))
		
		if (ncol(X)==1){
			stop("It is not possible to run cn.mops on only ONE sample.\n")
		}
		if (length(IRanges::unique(strand(input))) >1){
			stop(paste("Different strands found in GRanges object. Please make",
							"read counts independent of strand."))
		}
		chr <- as.character(seqnames(input))
		start <- start(input)
		end <- end(input)
		rownames(X) <- paste(chr,start,end,sep="_")
		
		irAllRegions <- IRanges(start,end)
		grAllRegions <- GRanges(chr,irAllRegions,seqinfo=seqinfo(input))
		grAllRegions <- sortSeqlevels(grAllRegions)
		names(irAllRegions) <- NULL
		
	} else if (is.matrix(input)){
		if (nrow(input)> 1){
			inputType <- "DataMatrix"
			X <- input
			X <- matrix(as.numeric(X),nrow=nrow(X))
			colnames(X) <- colnames(input)	
			chr <- rep("undef",nrow(X))
			irAllRegions <- IRanges(start=1:nrow(X),end=1:nrow(X))
			grAllRegions <- GRanges(chr,irAllRegions)
			
			
		} else{
			inputType <- "DataMatrix"
			chr <- "undef"
			irAllRegions <- IRanges(start=1:nrow(X),end=1:nrow(X))
			grAllRegions <- GRanges(chr,irAllRegions)
			parallel <- 0
		}
	} else if (is.vector(input)) {
		inputType <- "DataMatrix"
		X <- matrix(input,nrow=1)
		X <- matrix(as.numeric(X),nrow=nrow(X))
		chr <- "undef"
		irAllRegions <- IRanges(start=1:nrow(X),end=1:nrow(X))
		grAllRegions <- GRanges(chr,irAllRegions)
		parallel <- 0
	}else{
		stop("GRanges object or read count matrix needed as input.")
	}
	
	if (!all(isUnique(grAllRegions))) stop(paste("Genomic Ranges must be",
						"unique. Check \"all(isUnique(input))\" and remove identical segments."))
	
	
	if (any(X<0) | any(!is.finite(X))){
		stop("All values must be greater or equal zero and finite.\n")
	}
	if (!is.numeric(I)) stop("\"I\" must be numeric.")
	if (!is.character(classes)) stop("\"classes\" must be character.")
	if (length(I)!=length(classes)){
		stop("I and classes must have same length!")
	}
	if (!("CN1" %in% classes)){stop("One element of classes must be CN1 .\n")}
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
	params <- list("haplocn.mops",I,
			classes,
			priorImpact,cyc,
			normType,normQu,
			upperThreshold,lowerThreshold,
			minWidth,segAlgorithm,minReadCount,"CN1",version,paste(...))
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
	
	if (m < 100){
		warning(paste("For this small number of segments normalization",
						"might not be applicable."))
	}
	
	if (norm==0){
		X.norm <- X
		cov <- rep(1,N)
	} else if (norm==1) {
		message("Normalizing...")
		X.norm <- normalizeChromosomes(X,chr=chr,normType=normType,qu=normQu,
				sizeFactor=sizeFactor,quSizeFactor=quSizeFactor)
		cov <- rep(1,N)
	} else if (norm==2) {
		X.viz <- normalizeChromosomes(X,chr=chr,normType=normType,qu=normQu,
				sizeFactor=sizeFactor,quSizeFactor=quSizeFactor)
		X.norm <- X
		# robust estimates for the different coverages
		cov <- apply(X.norm,2,function(x) {
					mean(x[which(x>quantile(x,0.05) & x < quantile(x,0.95))])
				})
		if (median(cov)==0)
			stop("Median of the coverages is zero!")
		cov <- cov/median(cov)
	} else {
		stop("\"norm\" must be 0,1 or 2.")
	}
	params$cov <- cov
	
	
	message("Starting local modeling, please be patient...  ")
	
	if (parallel > 0){
		cl <- parallel::makeCluster(as.integer(parallel),type="SOCK")
		parallel::clusterEvalQ(cl,"haplocn.mopsCE")
	} 
	
	res <- list()
	
	for (chrom in chrOrder){
		message(paste("Reference sequence: ",chrom))
		chrIdx <- chrDf[chrom,1]:chrDf[chrom,2]
		
		cov <- rep(1,N)
		#n,N, I set
		idxCN1 <- which(classes=="CN1")
		alphaInit <- rep(0.05,n) 
		alphaInit[idxCN1] <- 0.6
		alphaInit <- alphaInit/ sum(alphaInit)
		alphaPrior <- rep(0,n)
		alphaPrior[idxCN1] <- priorImpact
		
		
		if (parallel==0){
			resChr <-apply(X.norm[chrIdx, ,drop=FALSE],1,haplocn.mopsCE, I=I,
					classes=classes,cov=cov,cyc=cyc,N=N,n=n,idxCN1=idxCN1,
					alphaInit=alphaInit,alphaPrior=alphaPrior,
					minReadCount=minReadCount)
		} else {
			resChr <- parallel::parApply(cl,X.norm[chrIdx, ,drop=FALSE],1,haplocn.mopsCE, I=I,
					classes=classes,cov=cov,cyc=cyc,N=N,n=n,idxCN1=idxCN1,
					alphaInit=alphaInit,alphaPrior=alphaPrior,
					minReadCount=minReadCount)
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
	CN <- t(sapply(res,.subset2,3))
	rownames(CN) <- rownames(X)
	colnames(CN) <- colnames(X)
	sINI <- t(sapply(res,.subset2,4))
	rownames(sINI) <- rownames(X)
	colnames(sINI) <- colnames(X)
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
		
			colnames(segDf) <- c("chr","start","end","mean","median","sample","CN")
			segDf <- segDf[ ,c("chr","start","end","sample","median","mean","CN")]			
			segDf <- segDf[order(match(segDf$chr,chrOrder),match(segDf$sample,colnames(X)),segDf$start), ]
					
	
	
			
	
	
			callsS <- matrix(NA,nrow=m,ncol=N)
			colnames(callsS) <- colnames(X)
			
			for (chrom in chrOrder){
				chrIdx <- chrDf[chrom,1]:chrDf[chrom,2]
				segDfTmp <- subset(segDf,chr==chrom)
				callsS[chrIdx, ] <- 
						matrix(rep(segDfTmp$mean,segDfTmp$end-segDfTmp$start+1),
								ncol=N)
			}
			
			
			
			
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
			segDf <- segDf[ ,c("chr","start","end","sample","median","mean","CN")]
			
			segDfSubset <- segDf[which(segDf$mean >= upperThreshold
									| segDf$mean <= lowerThreshold), ]	
			segDfSubset <- 
					segDfSubset[which((segDfSubset$end-segDfSubset$start+1)
											>= minWidth), ]
		}
		
		
		
		
		if (nrow(segDfSubset)>0){
			
			# Assembly of result object
			r <- new("CNVDetectionResult")
			cnvrR <- GenomicRanges::reduce(GRanges(seqnames=segDfSubset$chr,
							IRanges(segDfSubset$start,segDfSubset$end),
							seqinfo=seqinfo(grAllRegions)))
			cnvrR <- sortSeqlevels(cnvrR)
			
			cnvrCN <- matrix(NA,ncol=N,nrow=length(cnvrR))
			
			colnames(cnvrCN) <- colnames(X) 
			
			sampleNames <- segDfSubset$sample
			
			if (inputType=="GRanges"){
				ir <- IRanges()
				irCNVR <- IRanges()
				for (chrom in chrOrder){
					chrIdx <- chrDf[chrom,1]:chrDf[chrom,2]
					inputChr <- input[chrIdx]
					segDfSubsetChr <- subset(segDfSubset,chr==chrom)
					cnvrRChr <- cnvrR[which(as.character(
											seqnames(cnvrR))==chrom)]
					if (nrow(segDfSubsetChr) >0){
						ir <- c(ir,IRanges(start(inputChr)[
												segDfSubsetChr$start],
										end(inputChr)[segDfSubsetChr$end]))
						irCNVR <- c(irCNVR,IRanges(start(inputChr)[
												start(cnvrRChr)],
										end(inputChr)[end(cnvrRChr)]))
					}
				}
			} else if (inputType=="DataMatrix"){
				ir <- IRanges(start=segDfSubset$start,end=segDfSubset$end)
				irCNVR <- IRanges(start=start(cnvrR),end=end(cnvrR))
			}
			
			
			rd <- GRanges(seqnames=segDfSubset$chr,ir,
					seqinfo=seqinfo(grAllRegions),
					"sampleName"=sampleNames,
					"median"=segDfSubset$median,"mean"=segDfSubset$mean,
					"CN"=segDfSubset$CN)
			rd<- sortSeqlevels(rd)
			
			
			cnvr <- GRanges(seqnames=seqnames(cnvrR),irCNVR,
					seqinfo=seqinfo(grAllRegions),
					CN=cnvrCN)
			cnvr <- sortSeqlevels(cnvr)
			
			
			if (norm==2){
				r@normalizedData    <- X.viz		
			} else {
				r@normalizedData    <- X.norm			
			}		
			r@localAssessments  <- sINI
			r@individualCall   	<- callsS
			r@iniCall        	<- INI
			r@cnvs				<- rd
			r@cnvr				<- cnvr
			
			
			
			if (inputType=="GRanges"){
				irS <- IRanges()
				for (chrom in chrOrder){
					chrIdx <- chrDf[chrom,1]:chrDf[chrom,2]
					inputChr <- input[chrIdx]
					segDfChr <- subset(segDf,chr==chrom)
					if (nrow(segDfChr) >0 ){
						irS <- c(irS,IRanges(start(inputChr)[segDfChr$start],
										end(inputChr)[segDfChr$end]))
					}
				}
				r@segmentation <- GRanges(seqnames=segDf$chr,
						irS, seqinfo=seqinfo(grAllRegions),
						"sampleName"=segDf$sample,"median"=segDf$median,
						"mean"=segDf$mean,"CN"=segDf$CN)
				r@segmentation <- sortSeqlevels(r@segmentation) 
				
			} else if (inputType=="DataMatrix"){
				r@segmentation <- GRanges(seqnames=segDf$chr,
						IRanges(segDf$start,segDf$end),
						seqinfo=seqinfo(grAllRegions),
						"sampleName"=segDf$sample,
						"median"=segDf$median,
						"mean"=segDf$mean,"CN"=segDf$CN)
				
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
					inputChr <- input[chrIdx]
					segDfChr <- subset(segDf,chr==chrom)
					if (nrow(segDfChr) >0 ){
						irS <- c(irS,IRanges(start(inputChr)[segDfChr$start],
										end(inputChr)[segDfChr$end]))
					}
				}
				r@segmentation 			<- GRanges(seqnames=segDf$chr,
						irS, seqinfo=seqinfo(grAllRegions),
						"sampleName"=segDf$sample,"median"=segDf$median,
						"mean"=segDf$mean,"CN"=segDf$CN)
				r@segmentation <- sortSeqlevels(r@segmentation) 
				
			} else if (inputType=="DataMatrix"){
				r@segmentation 			<- GRanges(seqnames=segDf$chr,
						IRanges(segDf$start,segDf$end),
						seqinfo=seqinfo(grAllRegions),
						"sampleName"=segDf$sample,"median"=segDf$median,
						"mean"=segDf$mean,"CN"=segDf$CN)
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
}
