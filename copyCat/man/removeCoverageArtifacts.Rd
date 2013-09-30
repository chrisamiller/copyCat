\name{removeCoverageArtifacts}
\alias{removeCoverageArtifacts}
\title{
  remove alterations with extremely high or low coverage
}
\description{  
  identifies alterations with coverage at least 10x higher or lower than
  the median depth, which are generally artifacts. Removes them.
}
\usage{
  removeCoverageArtifacts(segs, rdo, rdo2)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{segs}{
    a data frame of altered segments returned by cnSegments() or
    cnSegments.paired(), then filtered to contain only alterations by
    getAlts().
  }
  \item{rdo}{
    a readDepth object created with new("rdObject") and filled with read
    counts using the getReadDepth() command (optionally corrected first)
  }
  \item{rdo2}{
    a readDepth object created with new("rdObject") and filled with read
    counts using the getReadDepth() command (optionally corrected
    first). This argument is used for the second in a pair of samples
    (i.e. matched tumor and normal). In other cases, it can be left blank.
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
  #  alts = removeCoverageArtifacts(segs,rdo)
}
