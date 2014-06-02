\name{subtractBins}
\alias{subtractBins}
\title{
  subtractBins
}
\description{
  takes two rd objects, and subtracts the read depth in each bin from the
  second from the first (using the rd column)
}
\usage{
  subtractBins(rdo,rdo2,normalization=FALSE)
}
\arguments{
  \item{rdo}{
    a readDepth object
  }
  \item{rdo2}{
    a readDepth object
  }
  \item{normalization}{
    do normalization of the read counts
  }
}
\value{
  an rd object identical to the first rd object, except that the read
  depth in the rd column is the difference between the rd columns from the two
  objects.
}

\examples{
  ##create a rdObject, set bin params
  # x <- new("rdObject")

  ##fill the bins with reads
  # x <- getReadDepth(x)

  ##create a rdObject, set bin params
  # y <- new("rdObject")

  ##fill the bins with reads
  # y <- getReadDepth(y)

  ##sum the rd columns
  # x <- subtractBins(x,y)
}
