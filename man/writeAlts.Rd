\name{writeAlts}
\alias{writeAlts}
\title{
  writeAlts
}
\description{
  write out a tab-delimited table of CN alterations to
  the file "output/alts.dat"
}
\usage{
  writeAlts(segs, rdo, filename="alts.dat")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{segs}{
    a data frame that is the result of calling rd.cnSegments
  }
  \item{rdo}{
    a readDepth object created with new("rdObject")
  }
  \item{filename}{
    the filename to which the output will be written
  }
}
\examples{

  ##get regions of CN gain and loss
  #  segs <- rd.cnSegments(x)

  ##write out the table
  #  writeAlts(segs, rdo)
}
