\name{cnSegments.paired}
\alias{cnSegments.paired}
\title{
  Identifies regions of differential copy-number alteration between two paired samples
}
\description{  
  
  uses Circular Binary Segmentation (as implemented in the DNAcopy
  package) to segment the genome and returns regions that contain CN
  gains or losses relative to the reference
}
\usage{
  cnSegments.paired(rdo.ref, rdo.test, onlyAlts=FALSE, minWidth=3, alpha=0.01, undoSD=2, rmGaps=TRUE)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{rdo.ref}{
    a readDepth object created from the reference sample using
    new("rdObject") and filled with read counts using the getReadDepth()
    command (optionally corrected first).
  }
  \item{rdo.test}{
    same as rdo, but containing the test sample's reads
  }
  \item{onlyAlts}{
    whether the function should return just those segments that 
    exceed the threshold values for gain and loss The Default is FALSE,
    which returns all segments.
  }
  \item{minWidth}{
    the minimum number of consecutive aberrant bins required before making
    a segment call. Default is three. 
  }
  \item{alpha}{
    the significance level for the CBS to accept change-points. Default 0.01
  }
  \item{undoSD}{
    require this many SDs between segment means to keep a split
  }  
  \item{rmGaps}{
    Controls whether adjacent segments have their boundaries adjusted to remove
    gaps. If so, the new segment boundary is set to halfway between the existing
    boundary calls.
  }
}
\value{
  returns a dataframe with 5 columns:
  chr  start  stop  #probes  copynumber

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

}
