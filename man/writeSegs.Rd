\name{writeSegs}
\alias{writeSegs}
\title{
  writeSegs
}
\description{
  write out a tab-delimited table of CN segments
  the file "output/segs.dat"
}
\usage{
  writeSegs(segs, rdo, filename="segs.dat")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{segs}{
    a data frame that is the result of calling rd.cnSegments
  }
  \item{rdo}{
    read depth object
  }
  \item{filename}{
    the filename to which the output will be written
  }

}
\examples{

  ##get regions of CN gain and loss
  #  segs <- rd.cnSegments(x)

  ##write out the table
  #  writeSegs(segs)
}
