\name{removeCoverageArtifacts}
\alias{removeCoverageArtifacts}
\title{
  remove alterations with extremely high or low coverage in normal sample
}
\description{  
  Identifies alterations with coverage at least 10x higher or lower than
  the median depth in the normal sample, which are generally artifacts.
  Removes them
}
\usage{
  removeCoverageArtifacts(segs, rdo)
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
    counts using the getReadDepth() command (optionally corrected
    first). This should be the object containing reads from the normal sample
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

  #  alts = removeCoverageArtifacts(segs,rdo)
}
