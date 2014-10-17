##-----------------------------------------------------------
## parse a 10-col mpileup file to extract SNP positions, depth
## and allelic fraction.
parseSamtools10col <- function(file, minimumDepth, maximumDepth){
  ## read in 10-col pileup+ file, grab only the pieces we need.
  sites = read.delim(file,header=F,quote="")
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

  sites$vaf = (sapply(sites$pileup,reflength)/sites$depth)*100

  ##drop the pileup strings to free up some memory
  sites$pileup <- NULL
  gc()

  ##throw out sites with depth greater than 100 and less than 20
  sites = sites[sites$depth >= minimumDepth & sites$depth <= maximumDepth,]
  sites$depth <- NULL
  
  return(sites)
}


##-----------------------------------------------------------
## parse a VCF file, extract SNP positions, depth, and
## allelic fraction
parseSamtoolsVcf <- function(file, minimumDepth, maximumDepth){
  #only read the columns we need (chr, pos, info)
  #assumes DP=\d+ is present in the info field
  sites = read.delim(file,header=F,quote="",comment.char="#", sep="\t",colClasses=c("character","integer",rep("NULL",5),"character","NULL","NULL"))
    
  names(sites) = c("chr","st","info")
  sites$vafs = as.numeric(str_replace(str_extract(sites$info,"AF1=[0-9.]+"),"AF1=",""))*100
    
  ## this would allow us to use DP4 field, which is slightly more accurate than DP, 
  ## but at the cost of 50x slower exection. We'll use the fast way for now.
  ## depths = unname(sapply(str_replace(str_extract(sites[,3],"DP4=[0-9,]+"),"DP4=",""),
  ##                        function(x){sum(scan(text=x,,sep=",",quiet=TRUE))}))
  
  depths = as.numeric(str_replace(str_extract(sites$info,"DP=[0-9]+"),"DP=",""))

  #subset by depth
  sites=sites[(depths >= minimumDepth & depths <= maximumDepth),]
  
  sites$info <- NULL
    
  return(sites)
}



##-----------------------------------------------------------
## infer the format of the samtools file
## 
inferSamtoolsFormat <- function(sfile){
  a = scan(file=sfile, sep="\n", nlines=10, what="character",quiet=T)
  if(grepl("fileformat=VCF",a[1],fixed=T)){
    return("VCF");
  }
  if(sum(grepl("##INFO",a,fixed=T)) > 0){
    return("VCF");
  }
  if(length(strsplit(a[1],"\t")[[1]]) == 10){
    return("10colPileup")
  }
  return("unknown")
}



##-----------------------------------------------------------
## read in a samtools pileup file (or VCF), calculate the VAF of each
## SNP, then find regions where the VAF is reasonably stable
## at 2x, 3x, or 4x. Use these regions to calculate the depth
## that is equivalent to CN-neutral regions
##
cnNeutralDepthFromHetSites <- function(rdo, samtoolsFile, snpBinSize, peakWiggle=3, minimumDepth=20, maximumDepth=100, plot=FALSE, samtoolsFileFormat="10colPileup",verbose=FALSE){

  library('stringr')

  #make sure we have a valid format
  recognizedFormats = c("10colPileup","VCF");  
  if(!(samtoolsFileFormat %in% recognizedFormats)){
    print(paste("unrecognized samtools file format given:",samtoolsFileFormat))
    print("attempting to infer the format from the header of the file")
    samtoolsFileFormat = inferSamtoolsFormat(samtoolsFile);
    if(!(samtoolsFileFormat %in% recognizedFormats)){
      print("unable to infer samtools format. Please try again providing one of the following recognized formats")
      print(paste(recognizedFormats,sep="\n"))
      stop("unable to continue");
    }
    print(paste("inferred that samtools file format is:",samtoolsFileFormat))
  }
  
  #read in the sites
  sites = c();
  if(verbose){ print("reading in samtools file")};
  
  if(samtoolsFileFormat=="10colPileup"){
    sites = parseSamtools10col(samtoolsFile, minimumDepth, maximumDepth)
  } else if (samtoolsFileFormat=="VCF"){
    sites = parseSamtoolsVcf(samtoolsFile, minimumDepth, maximumDepth)
  } else {
    stop(paste("ERROR: undefined samtools file format:",samtoolsFileFormat))
  }



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
