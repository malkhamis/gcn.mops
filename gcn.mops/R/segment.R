# Copyright (C) 2011 Klambauer Guenter 
# <klambauer@bioinf.jku.at>

#' @title Fast segmentation of CNV calls.
#' 
#' @description Performs a fast segmentation algorithm based on the cyber t test
#' and the t statistics. This is a special version for log-ratios or I/NI calls
#' that are assumed to be centered around 0. For segmentation of data with 
#' different characteristics you can a) substract the mean/median/mode from
#' your data or b) use the more general version of this algorithm in the R 
#' Bioconductor package "fastseg".
#' 
#' @param x Values to be segmented.
#' @param alpha Real value between 0 and 1 is interpreted as the percentage of
#' total points that are considered as initial breakpoints. An integer greater 
#' than 1 is interpreted as number of initial breakpoints. Default = 0.05.
#' @param segMedianT Vector of length 2. Thresholds on the segment's median. 
#' Segments' medians above the first element are considered as gains and below
#' the second value as losses. If set to NULL the segmentation algorithm tries
#' to determine the thresholds itself. If set to 0 the gain and loss segments
#' are not merged. (Default = NULL).
#' @param minSeg Minimum length of segments. Default = 3.
#' @param eps Real value greater or equal zero. A breakpoint is only possible 
#' between to consecutive values of x that have a distance of at least "eps".
#' Default = 0.
#' @param delta Positive integer. A parameter to make the segmentation more 
#' efficient. If the statistics of a breakpoint lowers while extending the 
#' window, the algorithm extends the windows by "delta" more points until it 
#' stops. Default = 20.
#' @param maxInt The maximum length of a segment left of the breakpoint and
#' right of the breakpoint that is considered. Default = 40.
#' @param cyberWeight The "nu" parameter of the cyber t-test. Default = 50.
#' @examples 
#' x <- rnorm(n=500,sd=0.5)
#' x[150:200] <- rnorm(n=51,mean=3,sd=0.5)
#' segment(x)
#' @author Guenter Klambauer \email{klambauer@@bioinf.jku.at}
#' @return A data frame containing the segments.
#' @importFrom IRanges as.data.frame
#' @importFrom IRanges setdiff
#' @importFrom IRanges width
#' @importFrom IRanges IRanges
#' @export
#' @useDynLib cn.mops


segment <- function(x, alpha=.05, segMedianT=NULL, minSeg=3, 
		eps=0, delta=20, maxInt=40, cyberWeight=50){
	
	if (any(!is.finite(x))){
		message("Detected infinite values in the data. Replacing with max/min!")
		y <- x[which(is.finite(x) & !is.na(x))]
		x[which(x==Inf)] <- max(y,na.rm=TRUE)
		x[which(x==-Inf)] <- min(y,na.rm=TRUE)
		
	}   
	#browser()
	
	#globalMedian <- median(x,na.rm=TRUE)
	#adjustment for log ratios or sI/NI calls
	globalMedian <- 0
	if (any(is.na(x))){
		message("NA values detected. Replacing with median.")
		x[is.na(x)] <- globalMedian 
	}
	
	if (is.null(segMedianT)) {
		segMedianT <- c()
		#segMedianT[1] <- mean(x, na.rm=TRUE)+2*sd(x, na.rm=TRUE)
		#segMedianT[2] <- mean(x, na.rm=TRUE)-2*sd(x, na.rm=TRUE)
		segMedianT[1] <- globalMedian+2*sd(x, na.rm=TRUE)
		segMedianT[2] <- globalMedian-2*sd(x, na.rm=TRUE)
		
	} else {
		if (length(segMedianT)==1){
			segMedianT <- c(abs(segMedianT), -abs(segMedianT))
		}
	}
	
	if (!is.numeric(alpha)) stop("\"alpha\" must be numeric!")
	if (!is.numeric(minSeg)) stop("\"minSeg\" must be numeric!")
	if (!is.numeric(maxInt)) stop("\"maxInt\" must be numeric!")
	if (!is.numeric(delta)) stop("\"minSeg\" must be numeric!")
	if (!is.numeric(cyberWeight)) stop("\"cyberWeight\" must be numeric!")
	
	if (minSeg < 2) minSeg <- 2
	if (maxInt < (minSeg+5)) maxInt <- minSeg+5
	if (cyberWeight < 0) cyberWeight <- 0

	m <- length(x)	
	
	if (missing("eps")){
		eps <- quantile(abs(diff(x)), probs=0.75)
	}
	
	res <-  .Call("segment", x, as.double(eps), as.integer(delta),
			as.integer(maxInt), as.integer(minSeg),
			as.double(cyberWeight))
	
	#message("Finished C function.")
	
	if (alpha >= 1){
		alpha <- as.integer(alpha)
		#message(paste("Number of initial breakpoints: ",alpha))
		brkptsInit <- sort(order(res$stat,decreasing=TRUE)[1:alpha])
		
	} else if (alpha < 1 & alpha > 0){
		#message(paste("Number of initial breakpoints: ", 
		#				as.integer(alpha*length(x)) ))
		pValT <- quantile(res$stat, probs=1-alpha)
		brkptsInit <- which(res$stat > pValT)
	} else{
		stop(paste("Alpha must be either between 0 and 1 or an integer",
						"greater than 1."))
	}
	
	nbrOfBrkpts <- length(brkptsInit)+1
	start <- c(1,brkptsInit+1)
	end <- c(brkptsInit,m)
	brkptsInit <- c(0, brkptsInit, m)
	
	
	avgs <- sapply(2:(nbrOfBrkpts+1),function(i){ 
						c(median(x[((brkptsInit[i-1]+1):brkptsInit[i])]),
								mean(x[((brkptsInit[i-1]+1):brkptsInit[i])]))
					}
			)
	
			
	df <- data.frame("start"=start, "end"=end, "mean"=avgs[2, ], 
			"median"=avgs[1, ])
	
	
	if (all(segMedianT==0)) {
		#message("No merging of segments.")
		ir <- IRanges::IRanges(df$start, df$end)
		ir <- ir[which(IRanges::width(ir)>=minSeg)]
		
		
		irAll <- IRanges::IRanges(1, length(x))
		segsFinal <- IRanges::as.data.frame(sort(c(ir, IRanges::setdiff(irAll, ir))))
		
		bIdx <- c(0,segsFinal$end)
		
		avgs <- sapply(2:(length(bIdx)),function(i){ 
					c(median(x[((bIdx[i-1]+1):bIdx[i])]),
							mean(x[((bIdx[i-1]+1):bIdx[i])]))
				}
		)
		
		df2 <- data.frame("start"=segsFinal$start, "end"=segsFinal$end, 
				"mean"=avgs[2,], "median"=avgs[1,])
		
		
		return(df2)
		
		
	} else {
		#browser()
		#message("Merging segments.")
		dfAmp <- df[which(df$median > segMedianT[1]), ]
		irAmp <- IRanges::IRanges(dfAmp$start, dfAmp$end)
		irAmp <- IRanges::reduce(irAmp)
		
		dfLoss <- df[which(df$median < segMedianT[2]), ]
		irLoss <- IRanges::IRanges(dfLoss$start, dfLoss$end)
		irLoss <- IRanges::reduce(irLoss)
		
		ir <- sort(c(irAmp, irLoss))
		ir <- ir[which(IRanges::width(ir)>=minSeg)]
		
		rm(irAmp, irLoss, dfAmp, dfLoss)    
		
		irAll <- IRanges(1, length(x))
		segsFinal <- as.data.frame(sort(
						c(ir, IRanges::setdiff(irAll, ir))))
		
		bIdx <- c(0,segsFinal$end)
		
		avgs <- sapply(2:(length(bIdx)),function(i){ 
					c(median(x[((bIdx[i-1]+1):bIdx[i])]),
							mean(x[((bIdx[i-1]+1):bIdx[i])]))
				}
		)
		
		df2 <- data.frame("start"=segsFinal$start, "end"=segsFinal$end, 
				"mean"=avgs[2, ], "median"=avgs[1, ])
		
		
		return(df2)
	}
}


