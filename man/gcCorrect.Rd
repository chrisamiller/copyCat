\name{gcCorrect}
\alias{gcCorrect}
\title{
  readDepth: GC-content correction
}
\description{
  takes a readDepth object containing read counts and uses loess
  correction to account for GC content bias. Updates the @chrs rd column
  and adds a "gc" column containing that bin's GC content percentage
}
\usage{
  gcCorrect(rdo, meth=FALSE, outlierPercentage=0.01, resolution = 0.001)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{rdo}{
    a readDepth object created with new("rdObject") and filled with read
    counts using the readDepth() command
  }
  \item{meth}{
    if this data consists of bisulfite reads and a methylation map is
    present in annotations/methMap, gc correction will be done while
    taking into account the bisulfite treatment that converts C -> U 
    and methylated bases that were protected from this conversion
  }
  \item{outlierPercentage}{
    Removes the highest and lowest N\% of reads before loess correction
    This prevents overfitting near the extremes, to which loess is 
    susceptible. Default value is 0.01
  }
  \item{resolution}{
    the resolution at which to bin the values (0.01=1\%)
  }
}
\details{
  Assumes that appropriate annotation files are in place (see full
  copyCat documentation for more information).
}
\value{
  returns a rd object with the read depth values altered to reflect
  gc content composition-induced bias.
}
\seealso{
  'getReadDepth'
}

\examples{

  ##create a rdObject, set bin params
  #  x <- new("rdObject")

  ##fill the bins with reads
  #  x.reads <- getReadDepth(x)

  ##do gc correction on the read counts
  #  x.gcCorr = gcCorrect(x)

}
