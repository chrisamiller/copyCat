##-----------------------------------------------------------
## read in a samtools pileup file, calculate the VAF of each
## SNP, then find regions where the VAF is reasonably stable
## at 2x, 3x, or 4x. Use these regions to calculate the depth
## that is equivalent to CN-neutral regions
##
cnNeutralDepthFromHetSites <- function(rdo, samtoolsFile, snpBinSize, peakWiggle=3, minimumDepth=20, maximumDepth=100, plot=FALSE){

  if(verbose){
    print("reading in samtools file")
    print(date())
  }

  ## read in 10-col pileup file, grab only the pieces we need.
  sites = read.delim(samtoolsFile,header=F,quote="")
  if(length(sites) < 10){
    print("expecting 10-column samtools pileup file - this file has less than 10 columns")
    stop();
  }
  sites = sites[,c(1,2,8,9)]
  names(sites) = c("chr","st","depth","pileup")
  sites$pileup = as.character(sites$pileup)
 
  ##return the number of bases that match the reference from a pileup string
  reflength <- function(s){
    return(unname(length(gregexpr('[,.]',s)[[1]])))
  }

  sites$ref = sapply(sites$pileup,reflength)
  sites$vaf = (sites$ref/sites$depth)*100

  ##drop the pileup strings to free up some memory
  sites$pileup <- NULL
  gc()

  ##throw out sites with depth greater than 100 and less than 20
  sites = sites[sites$depth >= minimumDepth & sites$depth <= maximumDepth,]

  ##throw out sites with vaf less than 15%
  ##too much noise around the margin there
  sites = sites[sites$vaf >= 15,]

  homsites = sites[which(sites$vaf >= 85),]
  hetsites = sites[which(sites$vaf < 85),]

  sites = NULL
  gc()

  if(verbose){
    print("finding clean windows of CN 2x, 3x, 4x")
  }

  doCalc <- function(rdo, snpBinSize, peakWiggle, chr, hetsites, homsites, plot){
    chrLen = getChrLength(chr,rdo@entrypoints)
    ## create an IRange for bins
    bins <- IRanges(start = (0:(ceiling(chrLen/snpBinSize)-1)*snpBinSize)+1, end = (1:ceiling(chrLen/snpBinSize))*snpBinSize)
    end(bins[length(bins)]) <- chrLen

    ##create IRange for dbsnp vafs
    homIsites = IRanges(start=homsites$st, end=homsites$st)
    hetIsites = IRanges(start=hetsites$st, end=hetsites$st)

    ##figure out which reads fall in which bins
    hetBinnedSites <- as.matrix(findOverlaps(bins,hetIsites))
    homBinnedSites <- as.matrix(findOverlaps(bins,homIsites))

    ##grab the VAFs, find the peaks, make the call
    findPeaks <- function(series,span=3){
      z <- embed(series, span)
      s <- span%/%2
      v<- max.col(z) == 1 + s
      result <- c(rep(FALSE,s),v)
      result <- result[1:(length(result)-s)]
      result
    }

    sts = c()
    sps = c()
    ploidy = c()
    reads = c()

    if(plot){
      dir.create(paste(rdo@params$outputDirectory,"/plots/vafplots",sep=""), showWarnings=FALSE)
      if(!(is.null(rdo@params$prefix))){
        pdf(file=paste(rdo@params$outputDirectory,"/plots/vafplots/",rdo@params$prefix,".vafs.",chr,".pdf",sep=""))
      } else {
        pdf(file=paste(rdo@params$outputDirectory,"/plots/vafplots/vafs.",chr,".pdf",sep=""))
      }
    }
    for(i in 1:length(bins)){

      #if we have sites to look at
      if(length(hetsites[which(hetBinnedSites[,1]==i),]$chr) > 0){

        ##exclude sites that look 1x (low het to homo ratio)
        hetcount = length(hetsites[which(hetBinnedSites[,1]==i),1])
        homcount = length(homsites[which(homBinnedSites[,1]==i),1])        
        if(hetcount/homcount < 0.75){
          cn = NULL
        } else {
          if(length(hetsites[which(hetBinnedSites[,1]==i),]$vaf) < 2){
            return(NA)
          } else {
            den=(density(hetsites[which(hetBinnedSites[,1]==i),]$vaf))
            peaks=den$x[findPeaks(den$y)]
            peakHeight=den$y[findPeaks(den$y)]
            
            ##filter out low, probably false peaks
            peaks = peaks[which(peakHeight > 0.01)]
            
            if(plot){
              plot(den, col="blue", main=paste(i," - ",peaks),xlim=c(0,100))
              abline(v=peaks)
              for(p in peaks){
                text(p+2,0,labels=round(p))
              }
            }

            ## remove low vaf noise and double peaks called by
            ## strangeness in peak calling function
            peaks = unique(sort(round(peaks[peaks>15])))
            cn = NULL
            
            if ((length(peaks) == 1) &
                (abs(peaks[1]-50) < peakWiggle)){
              cn = 2
              
            }else if ((length(peaks) == 2) &
                      (abs(peaks[1]-33.33) < peakWiggle) &
                      (abs(peaks[2]-66.66) < peakWiggle)){
              cn = 3
              
            } else if ((length(peaks) == 3) &
                       (abs(peaks[1]-25) < peakWiggle) &
                       (abs(peaks[2]-50) < peakWiggle) &
                       (abs(peaks[2]-75) < peakWiggle)){
              cn = 4
              
              ##often, we don't see the 50% peak because one of the
              ##two alleles is amplified twice
            } else if ((length(peaks) == 2) &
                       (abs(peaks[1]-25) < peakWiggle) &
                       (abs(peaks[2]-75) < peakWiggle)){
              cn = 4
            }
            
            ##store the valid sites for use
            if(!(is.null(cn))){
              
              if(plot){
                text(0,0,labels=paste("CN",cn))
              }
              
              sts = c(sts,(i-1)*snpBinSize)
              sps = c(sps,i*snpBinSize)
              ploidy = c(ploidy,cn)
            }
          }
        }
      }
    }
    if(plot){
      dev.off()
    }

    if(verbose){
      cat(paste("chr",chr," VAF anaylsis: \n  ",length(which(ploidy == 2)),"\t2x windows\n  ",length(which(ploidy == 3)),"\t3x windows\n  ",length(which(ploidy == 4)),"\t4x windows\n  ",length(bins)-length(ploidy),"\twindows ambiguous\n",sep=""))
    }

    ##now grab all the reads and match them up, adjust readcounts by ploidy, put the readcounts into a big list
    validWinds = IRanges(start=sts, end=sps)
    ##get all the positions and readcounts for this chr

    df = makeDf(rdo@chrs,rdo@params)
    df = df[which(df$chr==chr),]
    muts = IRanges(start=df$pos,end=df$pos+rdo@params$binSize)

    matches = as.matrix(findOverlaps(validWinds,muts))

    aReadDepths = c()
    ##take the readcounts in each valid window, adjust for ploidy and add to the list
    for(i in 1:length(validWinds)){
      reads = df[which(matches[,1]==i),3]

      ##adjust to diploid level, add to list
      aReadDepths = c(aReadDepths, median(reads / (ploidy[i]/2)))
    }
    if(length(reads) == 0){
      return(NA)
    }
    if(verbose){
      print(paste("mean:   ",mean(aReadDepths, na.rm=TRUE)))
      print(paste("median: ",median(aReadDepths, na.rm=TRUE)))
    }
    return(aReadDepths)
  }

  ## find locations of 1MB windows that are
  ## unambiguously 2x, 3x

  ##this uses a ton of memory, so we're just going to use one core for now
  options(cores = 1);
  mcoptions <- list(preschedule = FALSE)
  chr=NULL
  adjReadDepths = foreach(chr=rdo@entrypoints$chr, .combine="append",.options.multicore=mcoptions) %dopar% {
    doCalc(rdo, snpBinSize, peakWiggle, chr, hetsites[which(hetsites$chr==chr),], homsites[which(homsites$chr==chr),], plot=plot)
  }
  ##restore cores for future use
  options(cores = rdo@params$maxCores)

  
  if(length(adjReadDepths) <= 10){
    print("Too few interpretable windows were found in samtools data - falling back to use median read depth")
    return(getMedianReadCount(rdo))
  }
  
  ##return the average of the adj readcounts
  if(verbose){print(mean(adjReadDepths, na.rm=TRUE))}
  
  if(plot){
    pdf(file=paste(rdo@params$outputDirectory,"/plots/vafplots/means.pdf",sep=""))
    
    
    hist(adjReadDepths,breaks=100,col="darkgreen",xlim=c(0,(mean(adjReadDepths, na.rm=TRUE)*2)))
    mtext("red=mean,blue=med")
    abline(v=mean(adjReadDepths,na.rm=T),col="red")
    abline(v=median(adjReadDepths,na.rm=T),col="blue")
    dev.off()
  }
  
  if(verbose){print(date())};
  
  return(median(adjReadDepths, na.rm=TRUE));
}
