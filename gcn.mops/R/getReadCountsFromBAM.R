# Copyright (C) 2012 Klambauer Guenter 
# <klambauer@bioinf.jku.at>
.countBAM <- function(bamFile,sl,WL,mode,refSeqName,quiet=FALSE){		
	if (!quiet){message("Reading file: ",bamFile)}
	
	
	x <- c()
	if (length(bamFile)>1){
		stop("Single BAM file accepted as input.")
	}
	
	indexed <- TRUE		
	if (file.exists(paste(bamFile,".bai",sep=""))){
		bamIndex <- paste(bamFile,"",sep="")
	} else if (file.exists(gsub("bam$","bai",bamFile))) {
		bamIndex <- gsub(".bam$","",bamFile)
	} else {
		indexed <- FALSE
	}
	
	
	if (indexed){
		#message("Using the index of the BAM files.")
		for (i in 1:length(refSeqName)){
			refSeqTmp <- refSeqName[i]
			if (!quiet){message(paste(" ",refSeqTmp)) }
			gr <- GenomicRanges::GRanges(refSeqTmp,IRanges::IRanges(0,sl[i]))
			if (mode=="paired"){
				param <- Rsamtools::ScanBamParam(
						Rsamtools::scanBamFlag(isPaired = TRUE,
								isFirstMateRead=TRUE,
								isProperPair=TRUE),
						what=c("pos","mpos"),which=gr)
				#which=IRanges::RangesList(refSeqName))
				readPos <-Rsamtools::scanBam(bamFile,param=param,index=bamIndex)[[1]]
				pp <- ((readPos$pos+readPos$mpos)/2)
				gapSize <- abs(readPos$pos-readPos$mpos)
				rm("readPos")
				ppFiltered <- pp[!(gapSize>5*median(gapSize,na.rm=TRUE)|
									is.na(gapSize))]
				ppFiltered <- ppFiltered[!is.na(ppFiltered)]
			} else{
				param <- Rsamtools::ScanBamParam(Rsamtools::scanBamFlag(
								isPaired = FALSE),what=c("pos"),which=gr)
				#which=IRanges::RangesList(refSeqName))
				readPos <- Rsamtools::scanBam(bamFile,param=param,index=bamIndex)[[1]]
				ppFiltered <- readPos$pos[!is.na(readPos$pos)]
				rm("readPos")
			}
			if (length(ppFiltered)==0){
				warning("No reads found in file: ",bamFile,"\n",
						"You may have to change the mode to \"paired\" or",
						"\"unpaired\"")
			}
			if (sl[i]%%WL!=1){
				brkpts <- c(seq(1,sl[i],WL),sl[i])
			} else{
				brkpts <- seq(1,sl[i],WL)
			}
			if (any(ppFiltered>=sl[i])){
				warning(paste("Some read positions are greater than",
								"length of reference sequence! File: ",bamFile,"\n"))
				ppFiltered <- ppFiltered[ppFiltered<sl[i]]
			}
			
			x <- c(x,hist(ppFiltered,breaks=brkpts,plot=FALSE)$counts)
			
		}
		if (!quiet){message("\n") }
		
	} else {
		
		if (mode=="paired"){
			param <- Rsamtools::ScanBamParam(
					Rsamtools::scanBamFlag(isPaired = TRUE,
							isFirstMateRead=TRUE,
							isProperPair=TRUE),
					what=c("rname","pos","mpos"))
			#which=IRanges::RangesList(refSeqName))
			readPos <-Rsamtools::scanBam(bamFile,param=param)[[1]]
			
			for (i in 1:length(refSeqName)){
				refSeqTmp <- refSeqName[i]
				if (!quiet){message(paste(" ",refSeqTmp)) }
				
				readPosIdx <- (readPos$rname==refSeqTmp)
				tpos <- readPos$pos[readPosIdx]
				tmpos <- readPos$mpos[readPosIdx]
				pp <- ((tpos+tmpos)/2)
				gapSize <- abs(tpos-tmpos)
				ppFiltered <- pp[!(gapSize>5*median(gapSize,na.rm=TRUE)|
									is.na(gapSize))]
				ppFiltered <- ppFiltered[!is.na(ppFiltered)]
				
				if (length(ppFiltered)==0){
					warning("No reads found in file: ",bamFile,"\n",
							"You may have to change the mode to \"paired\" or",
							"\"unpaired\"")
				}
				if (sl[i]%%WL!=1){
					brkpts <- c(seq(1,sl[i],WL),sl[i])
				} else{
					brkpts <- seq(1,sl[i],WL)
				}
				if (any(ppFiltered>=sl[i])){
					warning(paste("Some read positions are greater than",
									"length of reference sequence! File: ",bamFile,"\n"))
					ppFiltered <- ppFiltered[ppFiltered<sl[i]]
				}
				
				x <- c(x,hist(ppFiltered,breaks=brkpts,plot=FALSE)$counts)
			}
			if (!quiet){message("\n")}
			
		} else{
			#message("No indices for the BAM files found.")
			
			param <- Rsamtools::ScanBamParam(Rsamtools::scanBamFlag(
							isPaired = FALSE),what=c("rname","pos"))
			#which=IRanges::RangesList(refSeqName))
			readPos <- Rsamtools::scanBam(bamFile,param=param)[[1]]
			
			
			for (i in 1:length(refSeqName)){
				refSeqTmp <- refSeqName[i]
				if (!quiet){message(paste(" ",refSeqTmp)) }
				
				
				readPosIdx <- (readPos$rname==refSeqTmp)
				tpos <- readPos$pos[readPosIdx] 
				ppFiltered <- tpos[!is.na(tpos)]
				
				if (length(ppFiltered)==0){
					warning("No reads found in file: ",bamFile,"\n",
							"You may have to change the mode to \"paired\" or",
							"\"unpaired\"")
				}
				if (sl[i]%%WL!=1){
					brkpts <- c(seq(1,sl[i],WL),sl[i])
				} else{
					brkpts <- seq(1,sl[i],WL)
				}
				if (any(ppFiltered>=sl[i])){
					#browser()
					warning(paste("Some read positions are greater than",
									"length of reference sequence! File: ",bamFile,"\n"))
					ppFiltered <- ppFiltered[ppFiltered<sl[i]]
				}
				
				x <- c(x,hist(ppFiltered,breaks=brkpts,plot=FALSE)$counts)
				
			}
			if (!quiet){message("\n")}
			
		}
	}
	
	return(x)
}


#' @title Calculation of read counts from BAM files.
#' @description Generates the read counts from BAM Files. 
#' These counts are necessary for CNV detection methods based
#' on depth of coverage information.
#' Note that the function is much faster, if the BAM files have an index file.
#' The index file is assumed to be in the same folder and have an identical
#' file name except that ".bai" is appended.
#' 
#' This function can also be run in a parallel version.
#' 
#' @param BAMFiles BAMFiles
#' @param sampleNames The corresponding sample names to the BAM Files. 
#' @param refSeqName Name of the reference sequence that should be analyzed.
#' The name must appear in the header of the BAM file. If it is not given
#' the function will select the first reference sequence that appears in the
#' header of the BAM files.
#' @param WL Windowlength. Length of the initial segmentation of the genome in
#' basepairs. Should be chosen such that on the average 100 reads are contained
#' in each segment. If not given, cn.mops will try to find an appropiate window 
#' length.
#' @param mode Possible values are "paired" and "unpaired", whether the mapping 
#' algorithm was using a "paired" or "unpaired" strategy.
#' @param parallel The number of parallel processes to be used for this function.
#' Default=0.
#' @examples 
#' BAMFiles <- list.files(system.file("extdata", package="cn.mops"),pattern=".bam$",
#' 	full.names=TRUE)
#' bamDataRanges <- getReadCountsFromBAM(BAMFiles,
#' 					sampleNames=paste("Sample",1:3),WL=5000,mode="unpaired")
#' X <- getReadCountsFromBAM(BAMFiles,
#' 					sampleNames=paste("Sample",1:3),WL=5000,mode="unpaired",parallel=2)
#' @return An instance of "GRanges", that contains the breakpoints of the 
#' initial segments and the raw read counts that were extracted from the BAM
#' files. This object can be used as input for cn.mops and other CNV detection
#' methods.
#' @importFrom Rsamtools scanBamHeader
#' @importFrom Rsamtools ScanBamParam
#' @importFrom Rsamtools scanBamFlag
#' @importFrom Rsamtools scanBam
#' @importFrom parallel makeCluster
#' @importFrom parallel clusterEvalQ
#' @importFrom parallel parApply
#' @importFrom parallel stopCluster 
#' @author Guenter Klambauer \email{klambauer@@bioinf.jku.at}
#' @export


getReadCountsFromBAM <- function(BAMFiles,sampleNames,refSeqName,WL,
		mode,parallel=0){
	
	if (missing(mode)){
		stop("Parameter \"mode\" must be \"paired\" or \"unpaired\"!")
	}
	
	if ((!(mode %in% c("paired","unpaired"))) ){
		stop("Mode parameter must be \"paired\" or \"unpaired\"!")
	}
	
	if (missing(sampleNames)){
		sampleNames <- basename(BAMFiles)	
	}
	
	headerInfo <- Rsamtools::scanBamHeader(BAMFiles)
	
	targets <- lapply(headerInfo,.subset2,1)
	sn <- names(targets[[1]])
	sl <- as.integer(targets[[1]])
	
	
	message(paste("Identified the following reference sequences: ",
					paste(unique(unlist(sn)),collapse=",")   ))
	
	if (missing(refSeqName)){
		refSeqName <- unique(unlist(sn))[1]
		message(paste("Missing \"refSeqName\"! Selecting",refSeqName,
						"as reference sequence."))
	} else{
		message("Using ",paste(refSeqName,collapse=", ")," as reference.")
	}
	
	if (all(file.exists(paste(BAMFiles,".bai",sep="")))){
		message("Using indexed BAM files.")
	} else if (all(file.exists(gsub("bam$","bai",BAMFiles)))){
		message("Using indexed BAM files.")
	} else {
		message("Note that this function is much faster, if the indices of ",
				"the BAM files are present.")
	}
	
	if (!(all(refSeqName %in% unique(unlist(sn))))){
		stop("RefSeqName does not match identified reference sequences.")
	}
	
	if (any(is.na(sn))){
		stop(paste(refSeqName,"does not appear in header!"))}
	
	if (length(targets)>1){
		for (i in 2:length(targets)){
			if (!(all(sn==names(targets[[i]])) & 
						all(sl==as.integer(targets[[i]]))  )){
				stop(paste("Targets in Header file of ",BAMFiles[i]," are not identical", 
								"to the header of the file",BAMFiles[1],"!"))
			} 
		}
	}
	
	
	nidx <- match(refSeqName,sn)
	sn <- sn[nidx]
	sl <- sl[nidx]
	
	sl <- as.integer(unique(sl))
	
	if (missing(WL)){
		message(paste("Missing \"WL\"! cn.mops will suggest an",
						"appropiate value for the window length."))
		sampleNames <- sampleNames[order(file.info(BAMFiles)$size)]
		BAMFiles <- BAMFiles[order(file.info(BAMFiles)$size)]
		xs <- sum(.countBAM(BAMFiles[1],sl=sl,WL=min(sl),mode=mode,
						refSeqName=refSeqName,quiet=TRUE))
		if (xs==0){
			message(paste("It is not possible calculate an appropiate",
							"window length - no reads map to that reference."))
			WL <- 25000
		} else{
			WL <- max(round(100*sum(sl)/(xs),-3),100)
		}
		message("Window length set to: ",WL)
	}
	
	if (parallel==0){
		XL <- lapply(BAMFiles,.countBAM,sl=sl,WL=WL,mode=mode,
				refSeqName=refSeqName)
		
	} else {
		message("Using parallel version of this function.")
		cl <- parallel::makeCluster(as.integer(parallel),type="SOCK")
		parallel::clusterEvalQ(cl,".countBAM")
		XL <- parallel::parLapply(cl,BAMFiles,.countBAM,sl=sl,WL=WL,mode=mode,
				refSeqName=refSeqName)	
		parallel::stopCluster(cl)
	}
	
	if (length(BAMFiles)==1){
		X <- as.matrix(unlist(XL),ncol=1)
	} else	{X <- do.call("cbind",XL)}
	colnames(X) <- BAMFiles
	
	ir <- IRanges::IRanges()
	chrv <- c()
	rn <- c()
	
	for (i in 1:length(refSeqName)){
		if (sl[i]%%WL!=1){
			brkpts <- c(seq(1,sl[i],WL),sl[i])
		} else{
			brkpts <- seq(1,sl[i],WL)
		}
#	browser()
		
		nSegm <- length(brkpts)
		#rn <- c(rn,paste(refSeqName[i],"_",brkpts[1:(nSegm-1)],"_",
		#				brkpts[2:(nSegm)]-1,sep=""))
		
		ir <- c(ir,IRanges::IRanges(start=brkpts[1:(nSegm-1)],
						end=brkpts[2:(nSegm)]-1))
		
		chrv <- c(chrv,rep(refSeqName[i],nSegm-1))
	}
	#rownames(X) <- rn
	
	#browser()
	mode(X) <- "integer"
	
	#gr <- GenomicRanges::GRanges(seqnames=chrv, ranges = ir,sampleNames=X)
	gr <- GenomicRanges::GRanges(seqnames=chrv, ranges = ir)
	
	colnames(X) <- sampleNames
	IRanges::values(gr) <- X
	#names(gr@elementMetadata@listData) <- sampleNames
	#IRanges::colnames(IRanges::elementMetadata(gr)) <- sampleNames
	
	
	gr <- sortSeqlevels(gr)
	GenomeInfoDb::seqlengths(gr) <- targets[[1]][match(seqlevels(gr),names(targets[[1]]))]
	
	return(gr)
}
