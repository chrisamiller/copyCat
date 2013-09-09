\name{mapCorrect}
\alias{mapCorrect}
\title{
  Correct read depth for mapability
}
\description{
  takes a readDepth object containing read counts and updates the rd
  counts, accounting for repetitive and unmappable genomic regions. Also
  adds a "map" column containing that bin's percentage mapability
}
\usage{
  mapCorrect(rdo, outlierPercentage=0.01, minMapability=0.60,
resolution=0.01, skipCorrection=FALSE)
}
\arguments{
  \item{rdo}{
    a readDepth object created with new("rdObject") and filled with read
    counts using the getReadDepth() command
  }
  \item{minMapability}{
    the minimum amount of mapability required. If the mapability of a
    region is below this level, the read depth is set to NA. 
  }
  \item{outlierPercentage}{
    Removes the highest and lowest N% of reads before loess correction
    This prevents overfitting near the extremes, to which loess is 
    susceptible. Default value is 0.01
  }
  \item{skipCorrection}{
    Set this to true to just remove low-mapability sites without
    doing loess correction on any of the readcounts. Helps prevent
    spurious hits in low-coverage regions (2 reads vs 1 read should
    not be interpreted as a CN doubling).
  }
  \item{resolution}{
    The amount of division between bins, expressed in percent. A value
    of 1 would bin regions with mapability percentage of
    [1,2,...,99,100].    
  }      
}
\details{
  Assumes that appropriate annotation files are in place.
}
\value{
  returns a readDepth object with the read depth values altered to
  remove bias due to mapability
}
\seealso{
  'getReadDepth'
}

\examples{
  ##create a rdObject, set bin params
  #  x <- new("rdObject")

  ##fill the bins with reads
  #  x <- getReadDepth(x)

  ##do mapping correction on the read counts
  #  x.mapCorr = mapCorrect(x)
}
