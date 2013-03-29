\name{cnNeutralDepthFromHetSites}
\alias{cnNeutralDepthFromHetSites}
\title{
  cnNeutralDepthFromHetSites
}
\description{
  examines heterozygous SNP positions in fixed windows to identify
  copy-number neutral regions. Finds the median read depth in these
  regions and uses it as the baseline for copy number estimation.  
}
\usage{
  cnNeutralDepthFromHetSites(rdo, samtoolsFile, snpBinSize, peakWiggle=3, minimumDepth=20, maximumDepth=100, plot=FALSE)
}
\arguments{
  \item{rdo}{
    a readDepth object filled with read counts
  }
  \item{samtoolsFile}{
    the path to a file containing the output of running 'samtools
    mpileup' on the bam file
  }
  \item{snpBinSize}{
    the size of the window to use for counting het sites and estimating
    peaks. Default of 1MB should be reasonable for most human genomes.
  }
  \item{peakWiggle}{
    The amount of "wiggle" allowed in classifying peaks.  For example,
    if peak wiggle is set to 5, a peak at 45% will be considered
    cn-neutral, while a peak at 44% will not.
  }
  \item{minimumDepth}{
    The minimum depth of coverage needed to consider a het snp
    site. Prevents sampling error due to low coverage.
  }
  \item{maximumDepth}{
    The maximum depth of coverage allowable at a het snp site. Prevents
    consideration of sites with aberrantly high depth due to mapping artifacts.
  }
  \item{plot}{
    Whether to generate density plots of each peak for visual
    review. Places these in the output directory, under plots/vafPlots
  }
}
\value{
  a number that represents the median depth of coverage in copy-number
  neutral sites.
}
\seealso{
}

\examples{
  # tumorSamtoolsFile = "samtools.mpileup.output" 
  # rdo@params$med = cnNeutralDepthFromHetSites(rdo, tumorSamtoolsFile,
  #                      snpBinSize=1000000, plot=TRUE)
}
