# Copyright (C) 2012 Klambauer Guenter 
# <klambauer@bioinf.jku.at>

#' @title Calculation of read counts from BAM files for predefined segments.
#' @description Generates the read counts from BAM Files for predefined segments. 
#' This is the appropiate choice for exome sequencing data, where the
#' bait regions, target regions or exons are the predefined segments.
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
#' @param GR A genomic ranges object that contains the genomic coordinates of
#' the segments. 
#' @param mode Possible values are "paired" and "unpaired", whether the mapping 
#' algorithm was using a "paired" or "unpaired" strategy. 
#' @param parallel The number of parallel processes to be used for this function.
#' Default=0.
#' @param BAIFiles The names of the BAI files that belong to the BAM files. The
#' vector has to be in the same order as the vector BAMFiles. If the BAI files have
#' the same name as the BAM files, only with ".bai" attached, this parameter needs
#' not be set. (Default = NULL).
#' @examples 
#' BAMFiles <- list.files(system.file("extdata", package="cn.mops"),pattern=".bam$",
#' 	full.names=TRUE)
#' gr <- GRanges(c("20","20"),IRanges(c(60000,70000),c(70000,80000)))
#' bamDataRanges <- getSegmentReadCountsFromBAM(BAMFiles,GR=gr,mode="unpaired")
#' bamDataRanges <- getSegmentReadCountsFromBAM(BAMFiles,GR=gr,mode="unpaired",parallel=2)
#' @return An instance of "GRanges", that contains the breakpoints of the 
#' initial segments and the raw read counts that were extracted from the BAM
#' files. This object can be used as input for cn.mops and other CNV detection
#' methods.
#' @importFrom Rsamtools countBam
#' @importFrom parallel makeCluster
#' @importFrom parallel clusterEvalQ
#' @importFrom parallel parApply
#' @importFrom parallel stopCluster 
#' @author Guenter Klambauer \email{klambauer@@bioinf.jku.at}
#' @export


getSegmentReadCountsFromBAM <- function(BAMFiles,GR,sampleNames,
		mode,parallel=0,BAIFiles=NULL){
	if (missing(mode)){
		stop("Parameter \"mode\" must be \"paired\" or \"unpaired\"!")
	}
	
	if ((!(mode %in% c("paired","unpaired"))) ){
		stop("Mode parameter must be \"paired\" or \"unpaired\"!")
	}
	
	if (missing(sampleNames)){
		sampleNames <- basename(BAMFiles)	
	}
	
	if (missing(GR) | !inherits(GR,"GRanges")){
		stop("You must submit the coordinates as GRanges object.")
	}
	
	
	if (!all(file.exists(paste(BAMFiles,".bai",sep=""))) & is.null(BAIFiles)){
		stop("The indices of the BAM files must be present.\n",
				"They are supposed to have the same file name with ",
				".bai appended.")
	} 
	
	message("This may take a couple of minutes per BAM file.",
			"Please be patient.\n\n")
	
	if (mode=="unpaired"){
		param <- Rsamtools::ScanBamParam(Rsamtools::scanBamFlag(
						isPaired = FALSE),which=GR)
	} else {
		param <- Rsamtools::ScanBamParam(Rsamtools::scanBamFlag(
						isPaired = TRUE),which=GR)
	}
	
	X <- matrix(NA,nrow=length(GR),ncol=length(BAMFiles))
	
	if (parallel==0){
		for (i in 1:length(BAMFiles)){
			message("Processing ",BAMFiles[i])
			if (is.null(BAIFiles)){
				X[,i] <- Rsamtools::countBam(BAMFiles[i],param=param)$records
			} else {
				newBAI <- gsub(".bai","",BAIFiles)
				X[,i] <- Rsamtools::countBam(BAMFiles[i],param=param,
						index=newBAI[i])$records				
			}
		}	
	} else {
		message("Using parallel version of this function.")
		cl <- parallel::makeCluster(as.integer(parallel),type="SOCK")
		parallel::clusterEvalQ(cl,"Rsamtools::countBam")
		if (is.null(BAIFiles)){
			XL <- parallel::parLapply(cl,BAMFiles,Rsamtools::countBam,param=param)
		} else {
			newBAI <- gsub(".bai","",BAIFiles)
			XL <- parallel::parLapply(cl,1:length(BAMFiles),function(jj){
						return(Rsamtools::countBam(BAMFiles[jj],param=param,index=newBAI[jj]))	
					})
						
		}
		
		parallel::stopCluster(cl)
		for (i in 1:length(BAMFiles)){
			X[,i] <- XL[[i]]$records
		}
		rm("XL")
	}
	
	
	
	colnames(X) <- BAMFiles
	
	mode(X) <- "integer"
	colnames(X) <- sampleNames
	
	IRanges::values(GR) <- X
#names(gr@elementMetadata@listData) <- sampleNames
#IRanges::colnames(IRanges::elementMetadata(GR)) <- sampleNames
#gr <- GRanges(GenomicRanges::seqnames(GR),IRanges::ranges(GR),
#		sampleNames=X)
	
	return(GR)
}
