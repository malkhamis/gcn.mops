### R code from vignette source 'cn.mops.Rnw'

###################################################
### code chunk number 1: cn.mops.Rnw:40-47
###################################################
options(width=75)
set.seed(0)
library(cn.mops)
library(Biobase)
library(GenomicRanges)
library(GenomeInfoDb)
cn.mopsVersion <- packageDescription("cn.mops")$Version


###################################################
### code chunk number 2: cn.mops.Rnw:140-141
###################################################
library(cn.mops)


###################################################
### code chunk number 3: cn.mops.Rnw:150-152 (eval = FALSE)
###################################################
## BAMFiles <- list.files(pattern=".bam$")
## bamDataRanges <- getReadCountsFromBAM(BAMFiles,mode="paired")


###################################################
### code chunk number 4: cn.mops.Rnw:156-157 (eval = FALSE)
###################################################
## res <- cn.mops(bamDataRanges)


###################################################
### code chunk number 5: cn.mops.Rnw:162-163 (eval = FALSE)
###################################################
## plot(res,which=1)


###################################################
### code chunk number 6: cn.mops.Rnw:167-172
###################################################
data(cn.mops)
resCNMOPS <- cn.mops(XRanges)
pdf("003.pdf")
plot(resCNMOPS,which=7,toFile=TRUE)
dev.off()


###################################################
### code chunk number 7: cn.mops.Rnw:232-236
###################################################
BAMFiles <- list.files(system.file("extdata", package="cn.mops"),pattern=".bam$",
		full.names=TRUE)
bamDataRanges <- getReadCountsFromBAM(BAMFiles,
		sampleNames=paste("Sample",1:3),mode="unpaired")


###################################################
### code chunk number 8: cn.mops.Rnw:241-242
###################################################
(bamDataRanges)


###################################################
### code chunk number 9: cn.mops.Rnw:261-263
###################################################
data(cn.mops)
ls()


###################################################
### code chunk number 10: cn.mops.Rnw:268-269
###################################################
head(XRanges[,1:3])


###################################################
### code chunk number 11: cn.mops.Rnw:272-273 (eval = FALSE)
###################################################
## resCNMOPS <- cn.mops(XRanges)


###################################################
### code chunk number 12: cn.mops.Rnw:277-278
###################################################
resCNMOPS <- calcIntegerCopyNumbers(resCNMOPS)


###################################################
### code chunk number 13: cn.mops.Rnw:284-285
###################################################
head(X[,1:3])


###################################################
### code chunk number 14: cn.mops.Rnw:288-289 (eval = FALSE)
###################################################
## resCNMOPSX <- cn.mops(X)


###################################################
### code chunk number 15: cn.mops.Rnw:293-294 (eval = FALSE)
###################################################
## resCNMOPSX <- calcIntegerCopyNumbers(resCNMOPSX)


###################################################
### code chunk number 16: cn.mops.Rnw:300-301 (eval = FALSE)
###################################################
## all(individualCall(resCNMOPSX)==individualCall(resCNMOPS))


###################################################
### code chunk number 17: cn.mops.Rnw:307-308 (eval = FALSE)
###################################################
## (resCNMOPS)


###################################################
### code chunk number 18: cn.mops.Rnw:312-313
###################################################
cnvs(resCNMOPS)[1:5]


###################################################
### code chunk number 19: cn.mops.Rnw:318-319
###################################################
cnvr(resCNMOPS)[1,1:5]


###################################################
### code chunk number 20: cn.mops.Rnw:325-326
###################################################
(CNVRanges[15,1:5])


###################################################
### code chunk number 21: cn.mops.Rnw:332-334
###################################################
ranges(cnvr(resCNMOPS))[1:2]
ranges(cnvr(resCNMOPS)) %over% ranges(CNVRanges)


###################################################
### code chunk number 22: cn.mops.Rnw:342-343 (eval = FALSE)
###################################################
## help(CNVDetectionResult)


###################################################
### code chunk number 23: cn.mops.Rnw:352-353 (eval = FALSE)
###################################################
## segplot(resCNMOPS,sampleIdx=13)


###################################################
### code chunk number 24: cn.mops.Rnw:356-359
###################################################
pdf("002.pdf")
segplot(resCNMOPS,sampleIdx=13,seqnames="chrA")
dev.off()


###################################################
### code chunk number 25: cn.mops.Rnw:378-379 (eval = FALSE)
###################################################
## plot(resCNMOPS,which=1)


###################################################
### code chunk number 26: cn.mops.Rnw:382-385
###################################################
pdf("001.pdf")
plot(resCNMOPS,which=1,toFile=TRUE)
dev.off()


###################################################
### code chunk number 27: cn.mops.Rnw:418-425 (eval = FALSE)
###################################################
## library(cn.mops)
## BAMFiles <- list.files(pattern=".bam$")
## segments <- read.table("targetRegions.bed",sep="\t",as.is=TRUE)
## gr <- GRanges(segments[,1],IRanges(segments[,2],segments[,3]))
## X <- getSegmentReadCountsFromBAM(BAMFiles,GR=gr,mode="unpaired")
## resCNMOPS <- exomecn.mops(X)
## resCNMOPS <- calcIntegerCopyNumbers(resCNMOPS)


###################################################
### code chunk number 28: cn.mops.Rnw:430-432
###################################################
resultExomeData <- exomecn.mops(exomeCounts)
resultExomeData  <- calcIntegerCopyNumbers(resultExomeData )


###################################################
### code chunk number 29: cn.mops.Rnw:435-436 (eval = FALSE)
###################################################
## plot(resultExomeData,which=5)


###################################################
### code chunk number 30: cn.mops.Rnw:439-442
###################################################
pdf("004.pdf")
plot(resultExomeData,which=5,toFile=TRUE)
dev.off()


###################################################
### code chunk number 31: cn.mops.Rnw:458-461 (eval = FALSE)
###################################################
## #the following removes the "chr" from reference sequence names
## library(GenomeInfoDb)
## seqlevels(gr) <- gsub("chr","",seqlevels(gr))


###################################################
### code chunk number 32: cn.mops.Rnw:467-469 (eval = FALSE)
###################################################
## gr <- GRanges(segments[,1],IRanges(segments[,2]-30,segments[,3]+30))
## gr <- reduce(gr)


###################################################
### code chunk number 33: cn.mops.Rnw:482-489
###################################################
resRef <- referencecn.mops(cases=X[,1],controls=rowMeans(X),
		classes=c("CN0", "CN1", "CN2", "CN3", "CN4", "CN5", "CN6",
				"CN7","CN8","CN16","CN32","CN64","CN128"),
		I = c(0.025, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 8, 16, 32, 64),
		segAlgorithm="DNAcopy")
resRef <- calcIntegerCopyNumbers(resRef)
(cnvs(resRef))


###################################################
### code chunk number 34: cn.mops.Rnw:506-508
###################################################
XchrX <- normalizeChromosomes(X[1:500, ],ploidy=c(rep(1,10),rep(2,30)))
cnvr(calcIntegerCopyNumbers(cn.mops(XchrX,norm=FALSE)))


###################################################
### code chunk number 35: cn.mops.Rnw:522-524 (eval = FALSE)
###################################################
## resHaplo <- haplocn.mops(X)
## resHaplo <- calcIntegerCopyNumbers(resHaplo)


###################################################
### code chunk number 36: cn.mops.Rnw:606-611 (eval = FALSE)
###################################################
## library(cn.mops); data(cn.mops)
## result <- calcIntegerCopyNumbers(cn.mops(XRanges))
## segm <- as.data.frame(segmentation(result))
## CNVs <- as.data.frame(cnvs(result))
## CNVRegions <- as.data.frame(cnvr(result))


###################################################
### code chunk number 37: cn.mops.Rnw:618-621 (eval = FALSE)
###################################################
## write.csv(segm,file="segmentation.csv")
## write.csv(CNVs,file="cnvs.csv")
## write.csv(CNVRegions,file="cnvr.csv")


###################################################
### code chunk number 38: cn.mops.Rnw:637-638 (eval = FALSE)
###################################################
## toBibtex(citation("cn.mops"))


