\name{removeGapSpanningSegments}
\alias{removeGapSpanningSegments}
\title{
  removes segments that span a reference assembly gap
}
\description{  
  intersects alterations with a bed file of reference assembly gap and
  removes regions that overlap, by default, at least 75\% with one of
  these gaps. This has proven to be a reliable indicator of
  false-positive calls.
}
\usage{
  removeGapSpanningSegments(segs,rdo,maxOverlap=0.75)
}

\arguments{
  \item{rdo}{
    a readDepth object created with new("rdObject") and filled with read
    counts using the getReadDepth() command (optionally corrected first)
  }
  \item{segs}{
    a data frame of altered segments returned by cnSegments() or
    cnSegments.paired(), then filtered to contain only alterations by
    getAlts().
  }
  \item{maxOverlap}{
    if a segment overlaps at least this much with a gap, it will be removed.
  }
}
\value{
  returns a dataframe with 5 columns:
  chr  start  stop  num.probes  copynumber

}
\seealso{
  'readDepth'
}
\examples{

  #create a rdObject, set bin params
  #  x <- new("rdObject")

  #fill the bins with reads
  #  x.reads <- readDepth(x)

  #do gc correction on the read counts
  #  x.gcCorr <- rd.gcCorrect(x)

  #return regions of CN gain and loss
  #  segs <- rd.cnSegments(x.gcCorr)

  #get alterations
  #  alts = getAlts(segs,rdo)

  #filter them
  #  alts = removeGapSpanningSegments(segs,rdo,maxOverlap=0.75)
}
