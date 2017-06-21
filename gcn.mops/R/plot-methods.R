# Copyright (C) 2011 Klambauer Guenter 
# <klambauer@bioinf.jku.at>

.convertToFastSegRes <- function(mopsres,segStat){
	ol <- IRanges::findOverlaps(gr(mopsres),segmentation(mopsres))
	ID <- as.character(values(segmentation(mopsres))$sampleName)
	if (segStat=="mean")
		seg.mean <- values(segmentation(mopsres))$mean
	else 
		seg.mean <- values(segmentation(mopsres))$median
	T <- table(subjectHits(ol))
	nn <- length(gr(mopsres))
	cs <- cumsum(T)%%nn
	num.mark <- as.integer(T)
	
	startRow <- c(1,(cs+1)[-length(cs)])
	endRow <- cs
	endRow[which(endRow==0)] <- nn
	grRet <- segmentation(mopsres)
	values(grRet) <- data.frame(ID=ID,num.mark=num.mark,seg.mean=seg.mean,
			startRow=startRow,endRow=endRow,stringsAsFactors=FALSE)
	
	return(grRet)
}

.makeLogRatios <- function(mopsres,mainCN="CN2"){
	X <- mopsres@normalizedData
	
	if (!is.null(mopsres@params$L[,mainCN])){
		r <- mopsres@params$L[,mainCN]
		#cat("lambda.\n")
	} else {
		#r <- apply(X,1,median)
		r <- Biobase::rowMedians(X)
	}
	
	R <- (X+0.01)/(r+0.01)
	#R[which(is.na(R))] <- 1
	
	#R[which(R==Inf)] <- max(R[which(is.finite(R))],na.rm=TRUE)
	
	#R[which(R==0)] <- min(R[which(R>0)],na.rm=TRUE)
	R <- log2(R)
	gr <- normalizedData(mopsres)
	colnames(R) <-  unique(as.character(
					IRanges::values(segmentation(mopsres))$sampleName))
	values(gr) <- R
	return(gr)
}


# A function to replace the sample names in a CNVDetectionResult
# oldNames <- sample(colnames(r@normalizedData))
# newNames <- sapply(strsplit(oldNames,"/|\\."),.subset2,8)
.replaceNames <- function(r, oldNames, newNames){
	#r is a CNVDetectionResult
	idx <- match(colnames(r@normalizedData),oldNames)
	colnames(r@normalizedData) <- newNames[idx]
	colnames(r@localAssessments) <-	newNames[idx]
	colnames(r@individualCall) <- newNames[idx]
	values(r@cnvs)["sampleName"] <-	newNames[match(as.character(values(r@cnvs)[,1]),oldNames)]
	CN <- values(r@cnvr)
	colnames(CN) <- newNames[idx]
	values(r@cnvr) <- CN
	values(r@segmentation)["sampleName"] <- 
			newNames[match(as.character(values(r@segmentation)[,1]),oldNames)]
	colnames(r@integerCopyNumber) <- newNames[idx]
	r@sampleNames <- newNames[idx]
	return(r)
}


#' @title Visualization of a CNV detection result.
#' 
#' @description Plots read counts, call values and CNV calls in an identified CNV region.
#' 
#' @param x An instance of "CNVDetectionResult" 
#' @param which The index of the CNV region to be plotted.
#' @param margin Vector of two positive integers that states how many segments 
#' left and right of the CNV region should be included in the plot. Default = 
#' c(10,10).
#' @param toFile Logical value whether the output should be plotted to a file.
#' Default = FALSE.
#' @examples
#' data(cn.mops)
#' r <- cn.mops(X[1:200, ])
#' plot(r)
#' @return Generates a CNV calling plot.
#' @author Guenter Klambauer \email{klambauer@@bioinf.jku.at}
#' @noRd
#' @export
#' @importFrom graphics plot
setMethod("plot", signature(x="CNVDetectionResult",y="missing"),
		function(x,which,margin=c(10,10),toFile=FALSE){
			
			if (missing(which)){
				message("Missing argument \"which\". Plotting the first CNVR.")
				which <- 1
			}
			
			tmp <- GRanges(seqnames(x@cnvr), IRanges(start(x@cnvr), start(x@cnvr)))
			MMstart <- findOverlaps(tmp, x@gr, select="first" )
			tmp <- GRanges(seqnames(x@cnvr), IRanges(end(x@cnvr), end(x@cnvr)))
			MMend <- findOverlaps(tmp, x@gr, select="first" )
			
			for (select in which){
				if (!toFile){
					dev.new()
				}
				
				if (select>length(x@cnvr)){
					stop("Selected unknown CNVR for plotting.")
				}
				
				refSeq <- as.character(seqnames(x@cnvr)[select])
				start <- MMstart[select]
				end <- MMend[select]
				plotStart <- max(start - margin[1], 1)
				plotEnd <- min(end + margin[2], nrow(x@normalizedData))
				
				if (plotEnd == plotStart) {
					plotStart <- plotStart - 1
				}
				
				layout(matrix(1:4,nrow=2))
				
				normDataSel <- x@normalizedData[plotStart:plotEnd, ,drop=FALSE]
				
				refSeqName <- unique(as.character(seqnames(x@cnvr)))
				xlab <- paste(refSeq,": ", unlist(start(x@gr))[plotStart], " - ",
						unlist(end(x@gr))[plotEnd], sep="")
				
				lt <- x@params$lowerThreshold
				ut <- x@params$upperThreshold
				col <- rep("grey", ncol(x@normalizedData))
				if (length(col)==1){
					sampOrd <- 1
					colOrd <- 1
				} else {
					if (is.numeric(lt) & is.numeric(ut)){
						col[apply(x@individualCall[start:end, ,drop=FALSE] >= ut, 2, any)] <- "red"
						col[apply(x@individualCall[start:end, ,drop=FALSE] <= lt, 2, any)] <- "blue"
					}
					sampOrd <- c(which(col == "grey"), which(col == "red"), which(col == "blue"))
					colOrd <- col[sampOrd]
				}
					
				
				lty <- sample(2:6, replace=TRUE, size=ncol(x@normalizedData))
				
				matplot(normDataSel[, sampOrd,drop=FALSE], type="l", lty=lty, lwd=2,
						main="Normalized Read Counts", ylab="Read Count",
						xlab=xlab, xaxt="n", col=colOrd)
				axis(1, at=1:length(plotStart:plotEnd),	labels=FALSE)
				
				matplot(x@localAssessments[plotStart:plotEnd, sampOrd,drop=FALSE],type="l",lty=lty,lwd=2,
						main="Local Assessments",ylab="Local Assessment Score",
						xlab=xlab,xaxt="n",col=colOrd)
				axis(1,at=1:length(plotStart:plotEnd), labels=FALSE)
				
				RA <- (normDataSel + 0.1) / apply(normDataSel + 0.1, 1, median)
				matplot(RA[, sampOrd,drop=FALSE],type="l",lty=lty,lwd=2,
						main="Read Count Ratios",ylab="Ratio",
						xlab=xlab,xaxt="n",col=colOrd)
				axis(1,at=1:length(plotStart:plotEnd),
						labels=FALSE)
				
				matplot(x@individualCall[plotStart:plotEnd, sampOrd,drop=FALSE],type="l",lty=lty,lwd=2,
						main="CNV Call",ylab="CNV Call Value",
						xlab=xlab,xaxt="n",col=colOrd)
				axis(1,at=1:length(plotStart:plotEnd),
						labels=FALSE)
			}
		})


#setMethod("segPlot", signature(r="CNVDetectionResult"),
setGeneric("segplot",
		function(r,mainCN="CN2", sampleIdx, seqnames, segStat="mean", plot.type="chrombysample", 
				altcol=TRUE, sbyc.layout, cbys.nchrom=1,
				cbys.layout, include.means=TRUE, zeroline=TRUE,
				pt.pch=".", pt.cex=3, pt.cols=c("green","black"),segcol, 
				zlcol="grey", ylim, lwd=3, ...) {
			standardGeneric("segplot")
		})

#' @title Visualization of a CNV detection result.
#' 
#' @description Plots the log normalized read counts and the detected segments as a 
#' segmentation plot.
#' 
#' @param r An instance of "CNVDetectionResult" 
#' @param mainCN The name of the main copy number. That is "CN2" for diploid
#' individuals. For haplocn.mops this should be set to "CN1".
#' @param sampleIdx The index of the samples to be plotted. (Default = missing)
#' @param seqnames The names of the reference sequence (chromosomes) to
#' be plotted. (Default = missing)
#' @param segStat Whether the segment line should display the mean or the 
#' median of a segments calls. (Default = "mean").
#' @param plot.type the type of plot. (Default = "s").
#' @param altcol logical flag to indicate if chromosomes should be
#'   plotted in alternating colors in the whole genome plot. (Default = TRUE).
#' @param sbyc.layout \code{layout} settings for the multifigure grid layout
#'    for the `samplebychrom' type.  It should be specified as a vector of
#'    two integers which are the number of rows and columns.  The default
#'    values are chosen based on the number of chromosomes to produce a
#'    near square graph.   For normal genome it is 4x6 (24 chromosomes)
#'    plotted by rows. (Default = NULL).
#' @param cbys.layout \code{layout} settings for the multifigure grid layout
#'   for the `chrombysample' type.  As above it should be specified as
#'    number of rows and columns and the default chosen based on the
#'    number of samples. (Default = NULL).
#' @param cbys.nchrom the number of chromosomes per page in the layout.
#'(Default = 1).
#' @param include.means logical flag to indicate whether segment means
#'   are to be drawn. (Default = TRUE).
#' @param zeroline logical flag to indicate whether a horizontal line at
#'    y=0 is to be drawn. (Default = TRUE).
#' @param pt.pch the plotting character used for plotting the log-ratio
#'    values. (Default = ".")
#' @param pt.cex the size of plotting character used for the log-ratio
#'    values (Default = 3).
#' @param pt.cols the color list for the points. The colors alternate
#'    between chromosomes. (Default = c("green","black").)
#' @param segcol the color of the lines indicating the segment means.
#' (Default = "red").
#' @param zlcol the color of the zeroline. (Default = "grey").
#' @param ylim this argument is present to override the default limits
#'    which is the range of symmetrized log-ratios. (Default = NULL).
#' @param lwd line weight of lines for segment mean and zeroline. (Default = 3).
#' @param ... other arguments which will be passed to \code{plot}
#'   commands.
#' @examples
#' data(cn.mops)
#' r <- cn.mops(X[1:200, ])
#' segplot(r,sampleIdx=1)
#' @return Generates a segmentation plot.
#' @author Guenter Klambauer \email{klambauer@@bioinf.jku.at}
#' @export


setMethod("segplot",
		signature(r="CNVDetectionResult"),
		function(r, mainCN="CN2",sampleIdx, seqnames, segStat="mean",
				plot.type="s", 
				altcol=TRUE, sbyc.layout, cbys.nchrom=1,
				cbys.layout, include.means=TRUE, zeroline=TRUE,
				pt.pch=".", pt.cex=3, pt.cols=c("green","black"),
				segcol="red", zlcol="grey", ylim, lwd=3, ...) {
			
			if (any(grepl("[^0-9A-Za-z]",colnames(r@normalizedData)))){
				message(
						paste("Segplot might not work because of special characters",
								"in the sample names. Use only A-Z,a-z and 0-9! \n There",
								"is a hidden function cn.mops:::.replaceNames that replaces",
								"the names in the \"CNVDetectionResult\" object."))
			}
			if (any(!grepl("^[A-Za-z]",colnames(r@normalizedData)))){
				message(
						paste("Segplot might not work because of special characters",
								"in the sample names. Use only A-Z,a-z and 0-9! \n There",
								"is a hidden function cn.mops:::.replaceNames that replaces",
								"the names in the \"CNVDetectionResult\" object."))
			}
			
			
			
			if (!missing(sampleIdx)){
				
				sn <- sampleNames(r)
				idx <- which(as.character(
								values(segmentation(r))
										$sampleName) %in% sn[sampleIdx])
				
				r@segmentation <- segmentation(r)[idx]
				#nd <- normalizedData(r)
				#IRanges::values(nd) <- IRanges::values(nd)[,sampleIdx]
				#IRanges::colnames(IRanges::values(nd)) <- sn[sampleIdx]
				r@normalizedData <- r@normalizedData[ ,sampleIdx,drop=FALSE]
				
				
			} 
			if (!missing(seqnames)){
				idx <- which(as.character(GenomicRanges::seqnames(
										segmentation(r))) %in% seqnames)
				r@segmentation <- segmentation(r)[idx]
				
				nd <- normalizedData(r)
				idx2 <- which(as.character(GenomicRanges::seqnames(nd))
								%in% seqnames)
				
				
				if (length(idx2)==0){
					stop(paste("Given \"seqnames\" do not appear in result",
									"object. Try to exchange \"chr1\" <--> \"1\"."))
				}
				#nd <- nd[idx2]
				
				r@gr <- r@gr[idx2]
				r@normalizedData <- r@normalizedData[idx2, ,drop=FALSE]
				if (!is.null(r@params$L)){
					r@params$L <- r@params$L[idx2, ]
				}
				
				
			} 
			R <- .makeLogRatios(r,mainCN)
			segm <- .convertToFastSegRes(r,segStat=segStat)			
			
			.segPlot(x=R,
					res=segm,
					plot.type=plot.type, 
					altcol=altcol, 
					cbys.nchrom=cbys.nchrom,
					include.means=include.means,zeroline=zeroline, 
					pt.pch=pt.pch,pt.cex=pt.cex, pt.cols=pt.cols,segcol=segcol,
					zlcol=zlcol, ylim=ylim, lwd=lwd, ...)
			
		})

