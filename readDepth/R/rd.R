##------------------------------------------------
## Given an rd object, add the number of reads in
## each bin.
readDepth <- function(rdo){

  ##set multicore options for reading if specified
  if("readCores" %in% names(rdo@params)){
    options(cores = as.numeric(rdo@params$readCores))
  }

  if(verbose){
    cat("Binning Start:",date(),"\n")
  }

  ##run bam-window here
    if (rdo@params$inputType == "bam") {
    rdo@params$binFile = paste(rdo@params$outputDirectory, "/bins.dat", sep="")
    getLibsFromBam(rdo@params)

    cmd = paste("/gscuser/cmiller/usr/src/bamwindow-v0.3/bam-window -q 1 -w ", rdo@params$inputBinSize, " -s ", sep="")
    if(rdo@params$perLib){
      cmd = paste(cmd,"-r ", rdo@params$outputDirectory, "/bamlibs ", sep="");
    }

    if (!(is.null(rdo@params$perLibLengths))){
      if (rdo@params$perLibLengths){
        cmd = paste(cmd, "-l ", sep="");
      }
    }
    cmd = paste(cmd, ">", rdo@params$binFile, sep="")
    
    if(verbose){
      print(cmd)
    }
    system(cmd)
  }
  

  if((rdo@params$inputType == "bins") | (rdo@params$inputType == "bam")) {
    ##read the bins in
    rdo = getWindowBins(rdo,rdo@params$binFile)

  } else { #bed files
    ##do the binning for each chr in parallel
    mcoptions <- list(preschedule = FALSE)
    b = foreach(chr=rdo@entrypoints$chr, .combine="combineBins",.options.multicore=mcoptions) %dopar% {
      doBinning(rdo@params, rdo@binParams, rdo@entrypoints, rdo@readInfo, chr)
    }
    rdo@chrs = b
  }
  closeAllConnections()
  gc() #just in case


  ## if we specified a bin size up front, go back and fill
  ## in the binParams information now (saves us from having
  ## read the file twice)
  if (!(is.na(rdo@params$binSize))){
    rdo@binParams = binParamsFromBins(rdo)
  }


  if(verbose){
    cat("Binning End:",date(),"\n")
    cat("-------------------------------\n")
    #plot distribution
    #plotDist(rdo, filename="output/dist.rawreads.pdf")
  }

  ##reset multicore options
  if("maxCores" %in% names(rdo@params)){
    options(cores = as.numeric(rdo@params$maxCores))
  }

  return(rdo)
}



##--------------------------------------------------
## do binning from a bin file (from bam-readcount output or from 
## pre-computed window counts)
getWindowBins <- function(rdo,binFile){
  winds = read.table(binFile,sep="\t",quote="", header=T)
  chrs = vector("list")

  #remove any columns with no data - all zeros
  for(i in names(winds)){
    if(sum(as.numeric(winds[[i]]))==0){
      winds[[i]] <- NULL
    }
  }
  
  ##if we're operating on a specific read length,
  ##grab just those columns
  if (!(is.null(rdo@params$perLibLengths))){
    if(rdo@params$perLibLengths){
      winds = winds[,c(1,2,grep(paste(".",rdo@params$readLength,"$",sep=""),names(winds)))]
    }
  }

  rdo@params$numLibs = length(winds)-2

  ##add the data to the rd object
  for(i in 1:length(rdo@entrypoints$chr)){
    chr = rdo@entrypoints$chr[i]
    
    tmp = data.frame(winds[which(winds[,1]==chr),1:length(winds)])
    names(tmp) = paste("rd.",names(tmp),sep="")

    ##handle the case where there's only one column, which tries to return a vector
    if(length(tmp)==3){
      rdo@chrs[[chr]] = data.frame(tmp[,3:length(tmp)])
      names(rdo@chrs[[chr]]) = names(tmp)[3]
    } else {
      rdo@chrs[[chr]] = tmp[,3:length(tmp)]
    }
  }

  return(rdo)
}


##--------------------------------------------------
## get the library and readgroup pairings from 
## the header of a bam file
getLibsFromBam <- function(params){
  library("digest")
  cat("creating bamlibs\n")
  header = scan(pipe(paste("samtools view -H",params$bamFile)), what="", sep="\n", quiet=TRUE)
  rglines = header[grepl("@RG",header)]  
  
  getRgLibPair <- function(line){
    a = strsplit(line,"\t")[[1]]
    rg = NULL;
    lib = NULL;

    if(grep("ID:",a[2])){
      rg = strsplit(a[2],":")[[1]][2];
    } else {
      return(NULL)
    }
    
    if(grep("LB:",a[5])){
      lib = strsplit(a[5],":")[[1]][2];
    } else {
      return(NULL)
    }

    libhash = digest(lib,algo="md5")
    
    return(paste(lib,rg,libhash,sep="\t"))
  }

  librg = lapply(rglines,getRgLibPair)
  
  write.table(librg,file=paste(params$outputDirectory,"/bamlibs",sep=""),sep="\n",quote=FALSE, row.names=FALSE, col.names=FALSE)

  
}

##--------------------------------------------------
## get the read depth for this chromosome and
## output a data frame in the proper position of the list
##
doBinning <- function(params, binParams, entrypoints, readInfo, chr){
  pref = paste(chr,":\t",sep="")

  if(verbose){
    cat(pref,"Starting\n")
  }

  binSize=binParams$binSize
  chrLen = getChrLength(chr,entrypoints)

  ## chunk the file if necessary to avoid running out of memory
  filename <- paste(chr,".bed",sep="")
  filepath <- paste(params$readDirectory,"/",chr,".bed",sep="")
  fileLen <- readInfo[which(readInfo$file==filepath),]$len

  b <- foreach(pos=seq(0,fileLen,by=params$chunkSize), .combine="+") %do% {
    binReads(binSize=binSize, chr=chr, chrLen=chrLen, filename=filepath, start=pos, chunkSize=params$chunkSize, zfileLen=fileLen)
  }
  if(verbose){
    cat(pref,"Done \n")
  }

  return(listify(data.frame(rd=b),chr))
}


##--------------------------------------------------
## map reads in bed format to bins for a given file
##
binReads <- function(binSize,chr,chrLen,filename,start,chunkSize,zfileLen){

  ## create an IRange for bins
  bins <- IRanges(start = (0:(ceiling(chrLen/binSize)-1)*binSize)+1, end = (1:ceiling(chrLen/binSize))*binSize)
  end(bins[length(bins)]) <- chrLen
  cat(chr,": reading chunk",(round(start/chunkSize)+1),"of",ceiling(zfileLen/chunkSize),"\n")

  ## input the reads
  rawreads = scan(pipe(paste("cut -f2",filename)), what=0, sep="\t", skip=start, nlines=chunkSize, quiet=TRUE)
  reads = IRanges(start = rawreads, end=rawreads)
  ##free up some memory
  fullreads = NULL

  ##figure out which reads fall in which bins
  binnedReads <- as.matrix(findOverlaps(bins,reads))

  nonZero = table(binnedReads[,1])

  a=rep(0,length(bins))
  b=names(nonZero)
  for(i in 1:length(b)){
    a[as.numeric(b[i])]=nonZero[i]
  }

  return(a)
}


##----------------------------------------------------
## create bin params from bins
##
binParamsFromBins <- function(rdo){
  bins = c()

  ##if this is a rdo with merged bins, there will be an "rd" column - use it:  
  foundrd = FALSE
  for(i in 1:length(names(rdo@chrs))){
    if(sum(grepl("^rd$",names(rdo@chrs[[names(rdo@chrs[i])]]))) > 0){
      bins = append(bins,rdo@chrs[[names(rdo@chrs)[i]]]$rd)
      foundrd=TRUE
    }
  }
  ##else, we need to add up the bins with "rd." prefixes indicating reads
  if(!(foundrd)){
    print("error: you must merge your libraries before calculating params - no changes made")
    return(rdo@binParams);
  }
  
  numReads = sum(bins)
  genomeSize = sum(rdo@entrypoints$length)
  mapPerc=sum(rdo@entrypoints$mapPerc * rdo@entrypoints$length) / sum(rdo@entrypoints$length)
  effectiveGenomeSize = genomeSize * mapPerc

  ploidyPerc = ploidyPercentages(effectiveGenomeSize,rdo@entrypoints, rdo@params)

  return(data.frame(binSize=rdo@params$binSize,
                    lossThresh = 1.5,  ##TODO - give params for these? use model to estimate?
                    gainThresh = 2.5,
                    med = (numReads/effectiveGenomeSize)*rdo@params$binSize,
                    hapPerc=ploidyPerc$haploidPerc,
                    dipPerc=ploidyPerc$diploidPerc,
                    tripPerc=ploidyPerc$triploidPerc))
}


##----------------------------------------------------
## add bins to the structure
## (from another library, for example)
##
addBins <- function(rdo,rdo2){
  for(chr in rdo@entrypoints$chr){
    if (length(rdo@chrs[[chr]]$rd) == length(rdo2@chrs[[chr]]$rd)){
      rdo@chrs[[chr]]$rd = rdo@chrs[[chr]]$rd + rdo2@chrs[[chr]]$rd
    } else {
      stop("bins don't match up properly")
    }
  }
  return(rdo)
}


##----------------------------------------------------
## subtract bins to the structure
## (from another library, for example)
##
subtractBins <- function(rdo,rdo2,normalization=FALSE){
  if(normalization){
    sum1 = 0
    sum2 = 0
    for(chr in rdo@entrypoints$chr){
      if (length(rdo@chrs[[chr]]$rd) == length(rdo2@chrs[[chr]]$rd)){
        sum1 = sum1 + sum(rdo@chrs[[chr]]$rd)
        sum2 = sum2 + sum(rdo2@chrs[[chr]]$rd)
        rdo@chrs[[chr]]$rd = rdo@chrs[[chr]]$rd - rdo2@chrs[[chr]]$rd
      } else {
        stop("bins don't match up properly")
      }
    }

    for(chr in rdo@entrypoints$chr){
      rdo@chrs[[chr]]$rd = rdo@chrs[[chr]]$rd + (rdo2@chrs[[chr]]$rd)*(sum1/sum2)
    }


  } else {
    for(chr in rdo@entrypoints$chr){
      if (length(rdo@chrs[[chr]]$rd) == length(rdo2@chrs[[chr]]$rd)){
        rdo@chrs[[chr]]$rd = rdo@chrs[[chr]]$rd + rdo2@chrs[[chr]]$rd
      } else {
        stop("bins don't match up properly")
      }
    }
  }

  return(rdo)
}

