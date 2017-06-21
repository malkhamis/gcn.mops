# Copyright (C) 2011 Klambauer Guenter 
# <klambauer@bioinf.jku.at>

# Modified by Mohammad Alkhamis <malkhamis@protonmail.com>
# this file is a modified version of 'cn.mops.R'

#cn.mops for external use for a single vector.
.cn.mopsC <- function(x,I = c(0.025,0.5,1,1.5,2,2.5,3,3.5,4), 
		classes=c("CN0","CN1","CN2","CN3","CN4","CN5","CN6","CN7","CN8"), cov,
		priorImpact = 1,cyc = 20,minReadCount=1) {
	
	N <- length(x)
	n <- length(I)
	
	if (missing(cov)){cov <- rep(1,N)}
	
	idxCN2 <- which(classes=="CN2")
	
	alpha.init <- rep(0.05,n)
	alpha.init[idxCN2] <- 0.6
	alpha.init <- alpha.init/ sum(alpha.init)
	alpha.est <- alpha.init
	alpha.prior <- rep(0,n)
	alpha.prior[idxCN2] <- 1
	alpha.prior <- alpha.prior*priorImpact
	
	if (all(x<=minReadCount)) {
		lambda.est <- rep(0,n)
		alpha.est <- rep(0,n)
		alpha.est[idxCN2] <- 1
		expCN <- rep(classes[idxCN2],N)
		ini <- 0
		ExpLogFoldChange <- rep(0,N)
		post.ik <- matrix(0,nrow=n,ncol=N)
		post.ik[idxCN2, ] <- 1
		
		
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
		ret=.Call("cnmops", as.numeric(x), I,cov, 
				as.integer(cyc), alpha.init, lambda.init,
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
				"params"=params)
		return(l)
	}
}


.segmentation <- function(x,chr,minWidth,DNAcopyBdry,...){
	m <- length(x)
	xx <- x+rnorm(mean=0,sd=0.00001,n=m)
	if (missing(chr))
		chr <- rep("undef",m)
	if (missing(minWidth))
		minWidth <- 4
	if (missing(DNAcopyBdry))	
		stop("\"DNAcopyBdry must be provided!")
	if (length(chr)!=m){
		stop("Vector \"chr\" must have the same length as \"x\"")
	}
	
	
	CNA.object <- DNAcopy::CNA(xx,
			chrom=chr,
			maploc=unlist(lapply(table(chr)[unique(chr)],function(x) 1:x)),
			data.type="logratio")
	
	
	segment.CNA.object <- DNAcopy::segment(CNA.object,
			min.width=minWidth, sbdry=DNAcopyBdry,...)
	segDf <- segment.CNA.object$output	
	names(segDf) <- c("sample","chr","start","end","idx","mean")
	#segDf$start <- pmax(segDf$start-1,1)
	#segDf$end <- c(segDf$end[1:(nrow(segDf)-1)]-1,segDf$end[nrow(segDf)])	
	
	#for calculating medians
	segDf$median <- apply(segDf,1,function(xx){
				median(x[which(chr==xx["chr"])][as.integer(xx["start"]):as.integer(xx["end"])])
			})
	return(segDf[,c("chr","start","end","mean","median")])
}


#' @title Copy number detection in NGS data. 
#' 
#' @description This function performs the cn.mops algorithm for copy number
#' detection in NGS data.
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
#' be assumed to have copy number 2. Default = 1.
#' @param cyc Positive integer that sets the number of cycles for the algorithm.
#' Usually after less than 15 cycles convergence is reached. Default = 20.
#' @param parallel How many cores are used for the computation. If set to zero
#' than no parallelization is applied. Default = 0.
#' @param gpu The GPU ID to be used for computation. If set to zero
#' then device ID zero is used. Default = 0.
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
#' applied to segments with very low coverage. Default=5.
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
#' cn.mops(XRanges)
#' cn.mops(XRanges,parallel=2)
#' 
#' @import methods
#' @import utils
#' @import stats
#' @importFrom parallel makeCluster
#' @importFrom parallel clusterEvalQ
#' @importFrom parallel parApply
#' @importFrom parallel stopCluster 
#' @importFrom GenomicRanges GRanges reduce values values<-
#' @importFrom IRanges as.matrix values unique IRanges reduce ranges
#' @importFrom S4Vectors subjectHits
#' @importFrom grDevices dev.cur dev.interactive dev.new
#' @importFrom GenomeInfoDb seqlevels sortSeqlevels seqnames seqinfo seqlengths
#' @importFrom BiocGenerics strand start end
#' @importFrom Biobase isUnique rowMedians
#' @importFrom graphics abline axis hist layout lines matplot mtext par title
# #' @importFrom stats density end median quantile rnorm sd start var
# #' @importFrom utils packageDescription
#' @useDynLib cn.mops
#' @return An instance of "CNVDetectionResult".
#' @author Guenter Klambauer \email{klambauer@@bioinf.jku.at}
#' @export

gcn.mops <- function(input,I = c(0.025,0.5,1,1.5,2,2.5,3,3.5,4),
		classes=c("CN0","CN1","CN2","CN3","CN4","CN5","CN6","CN7","CN8"),
		priorImpact = 1,cyc = 20,parallel=0, gpu=0,
		norm=1, normType="poisson",sizeFactor="mean",normQu=0.25, quSizeFactor=0.75,
		upperThreshold=0.5,lowerThreshold=-0.9,
		minWidth=3,segAlgorithm="fast",minReadCount=5,useMedian=FALSE,
		returnPosterior=FALSE,...){

############ check input ##################################################
# beginning time measurement
start.time <- Sys.time()

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

# beginning time measurement
start.time.chkX <- Sys.time()
	
	if (!all(isUnique(grAllRegions))) stop(paste("Genomic Ranges must be",
		"unique. Check \"all(isUnique(input))\" and remove identical segments."))
	
	if (any(X<0) | any(!is.finite(X))){
		stop("All values must be greater or equal zero and finite.\n")
	}
# end of time measurement
end.time.chkX <- Sys.time()
gpu.time.chkX <<- as.numeric(difftime(end.time.chkX ,start.time.chkX, unit = "secs"))

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
	
	# added by Alkhamis ##############
	if (!(is.numeric(gpu) & length(gpu)==1)) 
		stop("\"gpu\" must be numeric and of length 1.")
	##################################
	
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

# end of time measurement
end.time <- Sys.time()
gpu.time.formChkInput <<- as.numeric(difftime(end.time ,start.time, unit = "secs"))

	############################################################################
	
	version <- packageDescription("cn.mops")$Version
	params <- list("cn.mops",I,
			classes,
			priorImpact,cyc,
			normType,normQu,
			upperThreshold,lowerThreshold,
			minWidth,segAlgorithm,minReadCount,useMedian,"CN2",version,paste(...))
	names(params) <- c("method","folds",
			"classes",
			"priorImpact","cyc",
			"normType","normQu",
			"upperThreshold","lowerThreshold",
			"minWidth","segAlgorithm","minReadCount","useMedian","mainClass",
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
	
# beginning time measurement
start.time <- Sys.time()
	
	if (m < 100){
		warning(paste("Normalization might not be applicable",
						"for this small number of segments."))
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

# end of time measurement
end.time <- Sys.time()
gpu.time.normalization <<- as.numeric(difftime(end.time ,start.time, unit = "secs"))
	
# beginning of time measurement
start.time <- Sys.time()	

	message("Starting local modeling, please be patient...  ")

	message("Reference sequence(s) to process: ",paste(chrOrder, " "))

		idxCN2 <- which(classes=="CN2")
		alphaInit <- rep(0.05,n) 
		alphaInit[idxCN2] <- 0.6
		alphaInit <- alphaInit/ sum(alphaInit)
		alphaPrior <- rep(0,n)
		alphaPrior[idxCN2] <- priorImpact
		
		# Added/modified by alkhamis #####
		# we want to ensure that data are passed to the C function consistently. It must be a matrix 
		# (even if it has one row since R coerce a 1-row matrix into a vector. We are also converting 
		# the matrix to numeric (double). It was noted that if data were not normalized, R will coerce 
		# the matrix as "int" (GRange Object) while the C function expected "REAL()" is double.
		# This is a dirty workaround to satisfy "identical" funciton and this inconsistency should be corrected at the input level
		X.norm_numeric <- matrix(as.numeric(X.norm[, ,drop=FALSE]), ncol = ncol(X.norm[, ,drop=FALSE]))	
		res <- .Call("gcnmops_w", as.integer(gpu), X.norm_numeric,
		I, classes, cov, as.integer(cyc), idxCN2, alphaInit, alphaPrior, minReadCount,  returnPosterior)
		##################################
		
# end of time measurement
end.time <- Sys.time()
gpu.time.modeling <<- as.numeric(difftime(end.time, start.time, unit="secs"))
	
	#return(NULL)	
# beginning of time measurement
start.time <- Sys.time()
	
	## Postprocess result
	L <- res$lambda
	rownames(L) <- rownames(X)
	colnames(L) <- classes
	params$L <- L
	
	A <- res$alpha
	rownames(A) <- rownames(X)
	colnames(A) <- classes
	params$A <- A
	
	CN <- res$expectedCN
	rownames(CN) <- rownames(X)
	colnames(CN) <- colnames(X)
	
	sINI <- res$sini
	rownames(sINI) <- rownames(X)
	colnames(sINI) <- colnames(X)
	
	INI <- res$ini
	names(INI) <- rownames(X)
	
	# Modified by Alkhamis  ##########
	if (returnPosterior){
		post <- res$post
		tt <- try(attr(post, "dim") <- c(m,n,N))
		if (inherits(tt,"try-error")){
			message("Posterior too large for array extent.")
			post <- array(NA,dim=c(1,1,1))
		} else {
			dimnames(post) <- list(NULL,classes,colnames(X))
		}
	} else {
		post <- array(NA,dim=c(1,1,1))
	}
	##################################
	rm("res")
	
# end of time measurement
end.time <- Sys.time()
gpu.time.postprocessing <<- as.numeric(difftime(end.time, start.time, unit="secs"))

	#gpu_L <<- L
	#gpu_A <<- A
	#gpu_CN <<- CN
	#gpu_sINI <<- sINI
	#gpu_INI <<- INI
	#gpu_POST <<- post
	#warning("check point")
	#return(NULL)
	#res <- list("alpha_est"=A,"lambda_est"=L,"expCN"=CN,"sini"=sINI,"ini"=INI,"post"=post)
	#return(res)

	if (m>5){
start.time.sg <- Sys.time()
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
			
			#browser()
			
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
			
			if (useMedian){
				for (chrom in chrOrder){
					chrIdx <- chrDf[chrom,1]:chrDf[chrom,2]
					segDfTmp <- subset(segDf,chr==chrom)
					callsS[chrIdx, ] <- 
							matrix(rep(segDfTmp$median,segDfTmp$end-segDfTmp$start+1),
									ncol=N)
				}			
				segDfSubset <- segDf[which(
								segDf$median >= upperThreshold |
										segDf$median <= lowerThreshold), ]
			} else {
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
			}
		
			
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
				
				if (useMedian){
					callsS[chrIdx, ] <- 
							matrix(rep(segDfTmp$median,segDfTmp$end-segDfTmp$start+1),
									ncol=N)
				} else {
					callsS[chrIdx, ] <- 
							matrix(rep(segDfTmp$mean,segDfTmp$end-segDfTmp$start+1),
									ncol=N)
				}
	
				
				segDf <- rbind(segDf,segDfTmp)
			}
			
			#browser()
			segDf <- data.frame(segDf,"CN"=NA,stringsAsFactors=FALSE)		
			colnames(segDf) <- c("start","end","mean","median","sample",
					"chr","CN")
			segDf <- segDf[ ,c("chr","start","end","sample","median","mean","CN")]
			if (useMedian){
				segDfSubset <- segDf[which(segDf$median >= upperThreshold
										| segDf$median <= lowerThreshold), ]	
			} else {
				segDfSubset <- segDf[which(segDf$mean >= upperThreshold
										| segDf$mean <= lowerThreshold), ]	
			}
			segDfSubset <- 
					segDfSubset[which((segDfSubset$end-segDfSubset$start+1)
											>= minWidth), ]
		}
		
end.time.sg <- Sys.time()
gpu.time.seg <<- as.numeric(difftime(end.time.sg, start.time.sg, unit="secs"))

start.time <- Sys.time()		
		if (nrow(segDfSubset)>0){
			
			
			# Assembly of result object
			r <- new("CNVDetectionResult")
			cnvrR <- GenomicRanges::reduce(GRanges(seqnames=segDfSubset$chr,
							IRanges(segDfSubset$start,segDfSubset$end),
					seqinfo=seqinfo(grAllRegions)	))
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
			
			
			rd <- GRanges(seqnames=segDfSubset$chr,ir, seqinfo=seqinfo(grAllRegions),
					"sampleName"=sampleNames,
					"median"=segDfSubset$median,"mean"=segDfSubset$mean,
					"CN"=segDfSubset$CN)
			rd <- sortSeqlevels(rd)
			
			
			cnvr <- GRanges(seqnames=seqnames(cnvrR),irCNVR, seqinfo=seqinfo(grAllRegions))
			cnvr <-  sortSeqlevels(cnvr)
			GenomicRanges::values(cnvr) <- cnvrCN
			
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
				r@segmentation 			<- GRanges(seqnames=segDf$chr,
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
				r@segmentation <- sortSeqlevels(r@segmentation) 
				
			}
			
			r@gr <- grAllRegions
			r@posteriorProbs 	<- post
			r@params			<- params
			r@integerCopyNumber	<- CN
			r@sampleNames		<- colnames(X)

end.time <- Sys.time()
gpu.time.cnvCalling <<- as.numeric(difftime(end.time, start.time, unit="secs"))
			
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
						"mean"=segDf$mean,"CN"=segDf$CN )
				r@segmentation <- sortSeqlevels(r@segmentation) 
				
			} else if (inputType=="DataMatrix"){
				r@segmentation 			<- GRanges(seqnames=segDf$chr,
						IRanges(segDf$start,segDf$end), seqinfo=seqinfo(grAllRegions),
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

end.time <- Sys.time()
gpu.time.cnvCalling <<- as.numeric(difftime(end.time, start.time, unit="secs"))

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

end.time <- Sys.time()
gpu.time.cnvCalling <<- as.numeric(difftime(end.time, start.time, unit="secs"))
gpu.time.seg <<- 0

		return(r)	
		
	}
}
