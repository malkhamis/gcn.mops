\name{XRanges}
\docType{data}
\alias{XRanges}
\title{A simulated data set for CNV detection from NGS data.}
\description{This data set gives the read counts of 40 samples in 5000 genomic
locations. The rows correspond to genomic segments of 25kbp length and the
columns to samples. An entry is the number of reads that map to the specific
segment of the sample. The "GRanges" object contains the name of the 
reference sequence, start and end position of the genomic segments.
The simulated data contains CNVs given in the \code{\link{CNVRanges}}
object. It was generated using
distributions of read counts as they appear in real sequencing experiments. 
CNVs were implanted under the assumption that the expected read count is linear
dependent on the copy number (e.g. in a certain genomic we expect 
\eqn{\lambda}{lambda} reads for copy number 2, then we expect 
\eqn{2\cdot \lambda}{2*lambda} reads for copy number 4).
}
\usage{XRanges}
\format{A GRanges object with 5000 rows and 40 value columns across 1 space.}
\source{http://www.bioinf.jku.at/software/cnmops/cnmops.html.}
\references{Guenter Klambauer, Karin Schwarzbauer, Andreas Mayr, Djork-Arne Clevert, 
Andreas Mitterecker, Ulrich Bodenhofer, Sepp Hochreiter. 
\emph{cn.MOPS: mixture of Poissons for discovering copy number variations in next 
generation sequencing data with a low false discovery rate.}
Nucleic Acids Research 2012 40(9); doi:10.1093/nar/gks003.}
\keyword{datasets}
