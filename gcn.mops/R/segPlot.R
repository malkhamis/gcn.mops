# Plot method by Venkatraman E. Seshan, Adam Olshen (DNAcopy, v1.29)
# modified for fastseg input.


.segPlot <- function(x,res,
		plot.type, 
		altcol, sbyc.layout,
		cbys.nchrom, cbys.layout, include.means, zeroline,
		pt.pch, pt.cex, pt.cols, segcol, 
		zlcol, ylim, lwd,...){
	
	if (is.vector(x)){
		xdat <- as.matrix(x,ncol=1)
		nsample <- 1
		chrom <- as.character(seqnames(res))
		if (length(unique(chrom))!=1){
			stop("Data points and result indicate different number of chromosomes.")	
		}
		chrom <- rep(unique(chrom),length(x))
		maploc <- 1:length(x)
		xmaploc=FALSE 
	} else if (is.matrix(x)){
		if (is.null(colnames(x))){
			colnames(x) <- paste("Sample",1:ncol(x),sep="_")
		}
		xdat <- x
		nsample <- ncol(x)
		chrom <- as.character(seqnames(res))
		if (length(unique(chrom))!=1){
			stop("Data points and result indicate different number of chromosomes.")	
		}
		chrom <- rep(unique(chrom),nrow(x))
		maploc <- 1:nrow(x)
		xmaploc=FALSE 
	} else if (inherits(x, "GRanges")){
		#cat("..")	
		xdat <- do.call("cbind",values(x)@listData)
		nsample <- ncol(values(x))
		chrom <- as.character(seqnames(x))
		
		#chrom <- rep(unique(chrom),length(x))
		maploc <- (  start(ranges(x)) + end(ranges(x))  )/2
		xmaploc=TRUE
	}else{
		stop("x must be of type GRanges, vector or matrix!")
		
	}	
	
	if(missing(ylim)) {
		uylim <- max(abs(xdat[which(is.finite(xdat))]), na.rm=TRUE)
		ylim <- c(-uylim, uylim)
	}
	
	
	xres <- data.frame("ID" = values(res)$ID,
			"num.mark" = values(res)$num.mark,
			"chrom" = as.character(seqnames(res)),
			"seg.mean" = values(res)$seg.mean,stringsAsFactors=FALSE)
	xres <- xres[order(values(res)$ID,as.character(seqnames(res)),
					values(res)$startRow), ]
	
	#browser()
	
	
	
	if(dev.cur() <= 1) dev.new()
	int.dev <- dev.interactive()
	#plot.type <- match.arg(plot.type)
	op <- par(no.readonly = TRUE)
	parask <- par("ask")
	if (int.dev & !parask & nsample>1) par(ask = TRUE)
	sampleid <- colnames(xdat)
	if (is.null(sampleid)){sampleid <- paste("sample",1:nsample,sep="")}
	#browser()
	chrom0 <- chrom
	uchrom <- unique(chrom0)
	nchrom <- length(uchrom)
	if (xmaploc) {
		maploc0 <- as.numeric(maploc)
		if(length(uchrom)>1 & max(maploc0[chrom0==uchrom[1]]) > min(maploc0[chrom0==uchrom[2]])) {
			plen <- max(maploc0[chrom0==uchrom[1]])
			for(i in 2:nchrom) {
				#maploc0[chrom0==uchrom[i]] <- plen + maploc0[chrom0==uchrom[i]]
				maploc0[chrom0==uchrom[i]] <- maploc0[chrom0==uchrom[i]]
				plen <- max(maploc0[chrom0==uchrom[i]])
			}
		}
	}
	if (missing(pt.pch)) pt.pch <- "."
	if (missing(pt.cex)) {
		if (pt.pch==".") { pt.cex <- 3}
		else {pt.cex <- 1}
	}
	wcol0 <- rep(1, length(chrom0))
	if (altcol) {
		j <- 0
		for (i in uchrom) {
			j <- (j+1) %% 2
			wcol0[chrom0==i] <- 1+j
		}
	}
	
	if (missing(pt.cols)) pt.cols <- c("black","green")
	if (missing(segcol)) segcol <- "red"
	if (missing(zlcol)) zlcol <- "grey"
	if (missing(lwd)) lwd <- 3
	if (plot.type == "chrombysample" | plot.type == "c" ) {
		#cat("Setting multi-figure configuration\n")
		#browser()	
		par(mar = c(0, 4, 0, 2), oma = c(4, 0, 4, 0), mgp = c(2, 0.7, 0))
		if (missing(cbys.layout)) {
			nrow <- ncol <- ceiling(sqrt(nsample))
			if (nrow*ncol - nsample > 0) {
				nrow <- nrow - 1
				ncol <- ncol + 1
			}
			if (nrow*ncol - nsample >= nrow) ncol <- ncol - 1
			cbys.layout <- c(nrow, ncol)
		}
		
		lmat0 <- lmat1 <- c(1:nsample, rep(-cbys.nchrom*nsample,
						max(prod(cbys.layout) - nsample,1)))
		for(i in 1:(cbys.nchrom-1)) {
			lmat1 <- c(lmat1,lmat0+nsample*i)
		}
		lmat1[lmat1<0] <- 0
		lmat <- matrix(lmat1[1:prod(cbys.layout)], nrow = cbys.layout[1], 
				ncol = cbys.nchrom*cbys.layout[2], byrow = FALSE)
		#lmat <- matrix(lmat1, nrow = cbys.layout[1], ncol = cbys.layout[2], byrow = FALSE)
		#cat("lmat: ",lmat,"\n")
		#cat("cbys: ",cbys.layout,"\n")
		
		layout(lmat)
	}
	if (plot.type == "samplebychrom"| plot.type == "s" ) {
		#cat("Setting multi-figure configuration\n")
		par(mar = c(4, 4, 4, 2), oma = c(0, 0, 2, 0), mgp = c(2, 0.7, 0))
		if (missing(sbyc.layout)) {
			nrow <- ncol <- ceiling(sqrt(nchrom))
			if (nrow*ncol - nchrom > 0) {
				nrow <- nrow - 1
				ncol <- ncol + 1
			}
			if (nrow*ncol - nchrom > ncol) ncol <- ncol - 1
			sbyc.layout <- c(nrow, ncol)
		}
		lmat <- matrix(c(1:nchrom, rep(0,prod(sbyc.layout)-nchrom)),
				nrow = sbyc.layout[1], ncol = sbyc.layout[2], byrow=TRUE)
		layout(lmat)
	}
	if (plot.type == "chrombysample" | plot.type == "c") {
		atchrom <- 0.5/cbys.nchrom
		for (ichrom in uchrom) {
			if (xmaploc) maploc1 <- maploc0[chrom0==ichrom]
			for (isamp in 1:nsample) {
				genomdat <- xdat[chrom0==ichrom, isamp]
				ina <- which(is.finite(genomdat))
				genomdat <- genomdat[ina]
				if (xmaploc) maploc <- maploc1[ina]
				ii <- cumsum(c(0, 
								xres$num.mark[xres$ID == sampleid[isamp] & xres$chrom==ichrom]))
				mm <- xres$seg.mean[xres$ID == sampleid[isamp] & xres$chrom==ichrom]
				kk <- length(ii)
				zz <- cbind(ii[-kk] + 1, ii[-1])
				if (xmaploc) {
					plot(maploc, genomdat, pch = pt.pch, cex=pt.cex, 
							xaxt="n", ylim = ylim, ylab = sampleid[isamp])
				} else {
					plot(genomdat, pch = pt.pch, cex=pt.cex, xaxt="n",
							ylim = ylim, ylab = sampleid[isamp])
				}
				if(zeroline) abline(h=0, col=zlcol, lwd=lwd)
				if (isamp%%cbys.layout[1] == 0) {
					axis(1, outer=TRUE)
					title(xlab="Index")
				}
				if (include.means) {
					for (i in 1:(kk - 1)) {
						if (xmaploc) { 
							lines(maploc[zz[i, ]], rep(mm[i], 2), col = segcol, lwd=lwd)
						} else {
							lines(zz[i, ], rep(mm[i], 2), col = segcol, lwd=lwd)
						}
					}
				}
			}
			mtext(paste("Chromosome",ichrom), side = 3, line = 1, 
					at = atchrom, outer=TRUE, font=2)
			atchrom <- atchrom + 1/cbys.nchrom
			atchrom <- atchrom - floor(atchrom)
		}
	} else {
		for (isamp in 1:nsample)
		{
			#browser()
			xres <- xres[order(match(xres$chrom,uchrom)), ]
			genomdat <- xdat[, isamp]
			ina <- which(is.finite(genomdat))
			genomdat <- genomdat[ina]
			wcol <- wcol0[ina]
			chrom <- chrom0[ina]
			#if (xmaploc) maploc <- maploc0[ina]
			maploc <- maploc0[ina]
			ii <- cumsum(c(0, xres$num.mark[xres$ID == sampleid[isamp]]))
			mm <- xres$seg.mean[xres$ID == sampleid[isamp]]
			kk <- length(ii)
			zz <- cbind(ii[-kk] + 1, ii[-1])
			if(missing(ylim)) ylim <- range(c(genomdat, -genomdat))
			if (plot.type=="whole" | plot.type == "w"){
				#browser()
				#if (xmaploc) {
				#	plot(maploc, genomdat, pch = pt.pch, cex=pt.cex, 
				#			col=pt.cols[wcol], main = sampleid[isamp], 
				#			ylab = "", ylim = ylim)
				#	if(zeroline) abline(h=0, col=zlcol, lwd=lwd)
				#} else {
				plot(genomdat, pch = pt.pch, cex=pt.cex, 
						col=pt.cols[wcol], main = sampleid[isamp], 
						ylab = "", ylim = ylim)
				if(zeroline) abline(h=0, col=zlcol, lwd=lwd)
				#}
				if (include.means) {
					for (i in 1:(kk - 1))
					{
						#if (xmaploc) { 
						#	lines(maploc[zz[i, ]], rep(mm[i], 2),
						#			col = segcol, lwd=lwd)
						#} else {
						lines(zz[i, ], rep(mm[i], 2), col = segcol, lwd=lwd)
						#}
					}
				}
			}
			if (plot.type=="samplebychrom" | plot.type == "s")
			{
				#browser()
				
				cc <- xres$chrom[xres$ID == sampleid[isamp]]
				for (ichrom in uchrom)
				{
					if (xmaploc) {
						plot(maploc[chrom == ichrom], 
								genomdat[chrom == ichrom], pch = pt.pch, 
								cex=pt.cex, xlab="maploc", ylab = "", 
								main = paste("Chromosome", ichrom), ylim = ylim)
					} else {
						plot(genomdat[chrom == ichrom], pch = pt.pch, 
								cex=pt.cex, ylab = "", 
								main = paste("Chromosome", ichrom), ylim = ylim)
					}
					if(zeroline) abline(h=0, col=zlcol, lwd=lwd)
					if (include.means) {
						jj <- which(cc==ichrom)
						jj0 <- min(jj)
						for (i in jj)
						{
							if (xmaploc) {
								lines(maploc[zz[i, ]], rep(mm[i], 2), 
										col = segcol, lwd=lwd)
							} else {
								lines(1+zz[i, ]-zz[jj0,1], rep(mm[i], 2), 
										col = segcol, lwd=lwd)
							}
						}
					}
				}
				mtext(sampleid[isamp], side = 3, line = 0, outer = TRUE, font=2)
			}
			if (plot.type=="plateau" | plot.type == "p")
			{
				omm <- order(mm)
				ozz <- zz[omm,]
				ina <- unlist(apply(ozz, 1, function(ii) ii[1]:ii[2]))
				plot(genomdat[ina], pch = pt.pch, cex=pt.cex, 
						main = sampleid[isamp], ylab = "", ylim = ylim)
				if(zeroline) abline(h=0, col=zlcol, lwd=lwd)
				if (include.means) {
					ii <- cumsum(c(0, 
									xres$num.mark[xres$ID == sampleid[isamp]][omm]))
					smm <- mm[omm]
					zz <- cbind(ii[-kk] + 1, ii[-1])
					for (i in 1:(kk-1)) lines(zz[i, ], 
								rep(smm[i], 2), col = segcol, lwd=lwd)
				}
			}
		}
	}
}

