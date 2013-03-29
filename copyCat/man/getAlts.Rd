\name{getAlts}
\alias{getAlts}
\title{
  getAlts
}
\description{
  takes a data frame of segments and extracts just those containing
  copy-number gains and losses (that exceed the gain/loss thresholds)
}
\usage{
  getAlts(segs,rdo)
}
\arguments{
  \item{segs}{
    a 5-col data frame that is the result of calling cnSegments
  }
  \item{rdo}{
    a readDepth object
  }
}
\value{
  a data frame containing only segments of CN gain and loss
}
\seealso{
  'cnSegments'
}

\examples{

  ##get regions of CN gain and loss
  #  segs <- cnSegments(x)

  ##extract just the gains and losses
  #  alts <- getAlts(segs,x)
}
