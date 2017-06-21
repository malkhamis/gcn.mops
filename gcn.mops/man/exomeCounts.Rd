\name{exomeCounts}
\docType{data}
\alias{exomeCounts}
\title{Read counts from exome sequencing for CNV detection}
\description{This data set gives the read counts on chromosome 22 (hg19) 
of 22 samples in 3785 exons.  The rows correspond to targeted regions or exons 
and columns to samples.  An entry is the number of reads that map to the 
specific segment, i.e. targeted region or exon, of the sample. 
The GRanges object contains the information of the genomic location. The read 
counts were generated from freely available exome sequencing data of the 
1000Genomes Project.}
\usage{exomeCounts}
\format{A GRanges object of 3785 rows and 22 columns.}
\source{http://www.bioinf.jku.at/software/cnmops/cnmops.html.}
\references{Guenter Klambauer, Karin Schwarzbauer, Andreas Mayr, Djork-Arne Clevert, 
Andreas Mitterecker, Ulrich Bodenhofer, Sepp Hochreiter. 
\emph{cn.MOPS: mixture of Poissons for discovering copy number variations in next 
generation sequencing data with a low false discovery rate.}
Nucleic Acids Research 2012 40(9); doi:10.1093/nar/gks003.


The 1000 Genomes Project Consortium. 
\emph{A map of human genome variation from population-scale sequencing.}
Nature 2010 467(1061-1073); doi:10.1038/nature09534.}
\keyword{datasets}

