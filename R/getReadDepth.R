##------------------------------------------------
## Given an rd object, add the number of reads in
## each bin.
getReadDepth <- function(rdo){
  ##options to run bam-window from this package could go here

  rdo = getWindowBins(rdo)

  ## if(verbose){
  ##   cat("Binning End:",date(),"\n")
  ##   cat("-------------------------------\n")
  ##   ##plot distribution
  ##   ##plotDist(rdo, filename="output/dist.rawreads.pdf")
  ## }

  return(rdo)
}


##--------------------------------------------------
## create rdbins from a bin file (from bam-window output)
getWindowBins <- function(rdo){
  params=rdo@params
  winds = read.table(params$binFile,sep="\t",quote="", header=T)
  chrs = vector("list")

  rdo@params$binSize = winds[3,2]-winds[2,2]
  if(verbose){
    print(paste("inferred bin size: ", rdo@params$binSize))
  }

  ##remove any columns with no data - all zeros
  ## for(i in names(winds)){
  ##   if(sum(as.numeric(winds[[i]]))==0){
  ##     winds[[i]] <- NULL
  ##   }
  ## }

  ##sanity check
  if(length(winds) < 3){
    print("ERROR: the input file contains no bins with reads. Cannot continue")
    stop()
  }

  ##TODO - these checks should be improved
  ##if no per-lib information, can't do per-lib correction
  if((names(winds)[3] == "Counts") & (params$perLibrary == TRUE)){
    print("WARNING: bins file wasn't created with per-library option, will not perform per-library correction")
    params$perLibrary = FALSE;
  }
  ##ditto for per-read-length info
  if((names(winds)[3] == "Counts") & (params$perReadLength == TRUE)){
    print("WARNING: bins file wasn't created with per-read-length option, will not perform per-read-length correction")
    params$perReadLength = FALSE;
  }

  ## make sure we have a valid readlength from either the header or the params
  for(i in 3:length(winds)){
    splitname = strsplit(names(winds)[i],".",fixed=T)[[1]]
    if(is.na(suppressWarnings(as.integer(splitname[length(splitname)])))){
      if(params$readLength < 1){ #is set to zero as default
        print("Can't infer read length from window file headers")
        print("and read length was not specified in setParams()")
        print("correct one of these to proceed")
        stop()
      }
    }
  }


  ##create a list of all the libraries and read lengths
  createReadInfo <- function(winds, rdo){
    dataColNames = names(winds)[3:length(winds)]
    readInfo = NULL;
    for(i in 3:length(winds)){
      splitname = strsplit(names(winds)[i],".",fixed=T)[[1]]

      ##if readlength is specified in col name
      rl = c()
      if(!(is.na(suppressWarnings(as.integer(splitname[length(splitname)]))))){
        rl = as.integer(splitname[length(splitname)])
      } else { ##grab readlength from params
        rl = rdo@params$readLength
      }

      df = data.frame(column=i,
        lib=paste(splitname[1:(length(splitname)-1)],collapse="."),
        readlength=rl,
        numreads=sum(winds[,i]))

      if(is.null(readInfo)){
        readInfo <- df
      } else {
        readInfo <- rbind(readInfo,df )
      }
    }
    return(readInfo)
  }

  rdo@readInfo = createReadInfo(winds,rdo)


  ##we could have multiple libraries with the same read lengths
  ##if we aren't doing per-library correction, merge them
  mergeLibraries <- function(winds,rdo){
    torm=c()
    for(i in names(table(rdo@readInfo$readlength))){
      colsToMerge=rdo@readInfo[rdo@readInfo$readlength==i,]$column
      if(length(colsToMerge) > 1){
        newCol=rowSums(winds[,colsToMerge])
        winds=cbind(winds,newCol)
        names(winds)[length(winds)] = paste("Counts.",i,sep="")
        torm = c(torm, colsToMerge)
      }
    }
    if(length(torm)>1){
      for(j in rev(sort(torm))){
        winds[[j]] <- NULL
      }
    }
    return(winds)
  }

  ## we could have multiple read lengths from the same library
  ## if we aren't doing per-readlength correction, merge them
  mergeReadLengths <- function(winds,rdo){
    torm=c()
    for(i in names(table(rdo@readInfo$lib))){
      colsToMerge=rdo@readInfo[rdo@readInfo$lib==i,]$column
      if(length(colsToMerge) > 1){
        newCol=rowSums(winds[,colsToMerge])
        winds=cbind(winds,newCol)
        names(winds)[length(winds)] = paste(i,".lib",sep="")
        torm = c(torm, colsToMerge)
      }
    }
    if(length(torm)>1){
      for(j in rev(sort(torm))){
        winds[[j]] <- NULL
      }
    }
    return(winds)
  }



  ##if we're not doing per-library or per-readlength correction,
  ##just merge all of the reads appropriately
  if(!(params$perLibrary)) {
    winds = mergeLibraries(winds, rdo)
    rdo@readInfo = createReadInfo(winds,rdo)
  }
  if(!(params$perLibrary)) {
    winds = mergeReadLengths(winds,rdo)
    rdo@readInfo = createReadInfo(winds,rdo)
  }



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


  #infer bin size and calculate overall medan
  med = c()
  if(length(winds) > 3){
    med = median(rowSums(winds[,c(3:length(winds))]),na.rm=T)
  } else {
    med = median(winds[,3],na.rm=T)
  }

  rdo@binParams <- data.frame(binSize = winds[2,2]-winds[1,2],
                              med = med,
                              numLibs=length(names(tmp))-2
                              )

  return(rdo)
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
## subtract bins from the structure
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

