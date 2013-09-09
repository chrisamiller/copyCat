##-------------------------------------------------
## gc correction functions
##
gcCorrect <- function(rdo, meth=FALSE, outlierPercentage=0.01, resolution=0.001){
  if(verbose){
    cat("correcting for GC bias",date(),"\n")
  }

  ##have to correct each read length individually
  for(len in unique(rdo@readInfo$readlength)){
    if(verbose){
      print(paste("calculating GC content for readlength",len,date(),sep=" "))
    }
    
    rdo2 = subsetByReadLength(rdo,len)
    rdo2@params$annotationDirectory = getAnnoDir(rdo2@params$annotationDirectory, len)

    
    ##figure out the avg num of reads for each level of GC content
    mcoptions <- list(preschedule = FALSE)
    chr = NULL;
    gcBins = foreach(chr=rdo2@entrypoints$chr, .combine="combineBins",.options.multicore=mcoptions) %dopar% {
      binGC(rdo2, chr, len)
    }
    for(i in rdo2@entrypoints$chr){
      rdo2@chrs[[i]] = cbind(rdo2@chrs[[i]],gcBins[[i]])
    }
    closeAllConnections()

    ##if bisulfite reads, consider methylation before correcting
    if(meth){
      if(verbose){cat("Correcting GC content for bi-sulfite treatment and methylation:  ",date(),"\n") }
      foreach(chr=rdo2@entrypoints$chr) %dopar%{
        rdo2@chrs[[chr]] = methAdjust(rdo2,chr)
      }
    }

    ##finally, do the loess correction
    if(verbose){cat("Correcting read depth for GC-content bias:  ",date(),"\n")}

    ##for each library
    for(i in 1:(length(names(rdo2@chrs[[1]]))-1)){
      name = names(rdo2@chrs[[1]][i])
      if(verbose){
        print(paste("correcting library ",i," (",name,")",sep=""))
      }
      gcAdj = loessCorrect(rdo2,i,len,name,outlierPercentage=outlierPercentage,type="gc", corrResolution=resolution)
      for(j in rdo2@entrypoints$chr){
        rdo2@chrs[[j]][[name]]= gcAdj[[j]]
      }
    }

    ##  if(verbose){
    ##    plotDist(rdo2, filename="output/dist.postGcCor.pdf")
    ##    cat("Done",date(),"\n")
    ##    cat("------------------------\n")
    ##  }

    ##merge back into rdo
    rdo = replaceReadCounts(rdo,rdo2)
  }

  return(rdo)
}



##-------------------------------------------------
## mapability correction functions
##
mapCorrect <- function(rdo, outlierPercentage=0.01, minMapability=0.60, resolution=0.01, skipCorrection=FALSE){

  ##have to correct each read length individually
  for(len in unique(rdo@readInfo$readlength)){
    if(verbose){
      print(paste("calculating mapability content for readlength",len,date(),sep=" "))
    }

    rdo2 = subsetByReadLength(rdo,len)
    rdo2@params$annotationDirectory = getAnnoDir(rdo2@params$annotationDirectory, len)

    ##figure out the avg num of reads for each level of mapability
    mcoptions <- list(preschedule = FALSE)
    mapBins = foreach(chr=rdo2@entrypoints$chr, .combine="combineBins",.options.multicore=mcoptions) %dopar% {
      binMap(rdo2, chr, len)
    }
    for(i in rdo2@entrypoints$chr){
      rdo2@chrs[[i]] = cbind(rdo2@chrs[[i]],mapBins[[i]])
    }
    closeAllConnections()

    ##find bins with mapability below threshold, set value to NA
    for(chr in rdo2@entrypoints$chr){
      for(i in 1:rdo2@binParams$numLibs){
        if(length(which(rdo2@chrs[[chr]]$map < minMapability)) > 0){
          rdo2@chrs[[chr]][which(rdo2@chrs[[chr]]$map < minMapability),i] = NA
        }
      }
    }

    ##skipping correction can be used to just remove low-mapability sites. Helps prevent
    ##spurious hits in low-coverage regions. (2 reads vs 1 read should not be significant)
    if(skipCorrection==FALSE){
      ##finally, do the loess correction
      if(verbose){cat("Correcting read depth for mapability bias:  ",date(),"\n")}
      
      ##for each library
      for(i in 1:(length(names(rdo2@chrs[[1]]))-1)){
        if(verbose){cat("correcting library ",i,"\n");}
        name = names(rdo2@chrs[[1]][i])
        mapAdj = loessCorrect(rdo2,i,len,name,outlierPercentage=outlierPercentage, type="map", corrResolution=resolution)
        for(j in rdo2@entrypoints$chr){
          rdo2@chrs[[j]][[name]]= mapAdj[[j]]
        }
      }
    } 

    ##merge back into rdo
    rdo = replaceReadCounts(rdo,rdo2)
  }
  return(rdo)
}


##---------------------------------------------------
## adjust gc content to account for bisulfite sequencing
## and methylated bases (that don't get bs-converted)
methAdjust <- function(rdo,chr){
  readBin = rdo@chrs[[chr]]
  params = rdo@binParams
  entrypoints = rdo@entrypoints
  ##start by halving gc content (assume all C -> U)
  binSize = params$binSize
  if(!is.null(readBin$map)){
    effectiveLength = readBin$map * binSize
  }else{
    effectiveLength=binSize
  }
  gcBases <- readBin$gc/2 * effectiveLength

  ##now, figure out how many protected (methylated) bases were in each bin
  chrLen = getChrLength(chr,entrypoints)
  bins <- IRanges(start = (0:(ceiling(chrLen/binSize)-1)*binSize)+1, end = (1:ceiling(chrLen/binSize))*binSize)
  end(bins[length(bins)]) <- chrLen

  filename = paste(params$annotationDirectory,"/methMap/map.",chr,".gz",sep="")

  ## TODO - make this happen without calling awk
  ##use awk to quickly strip out lines that don't contain any methylated bases.
  ##this is more than an order of magnitude faster than doing it in R.
  system(paste("gunzip -c ",filename," | cut -f 1,4,6 | awk '$2 != \"0\"' >",filename,".tmp",sep=""))

  f = file(paste(filename,".tmp",sep=""))
  open(f)

  ##get the total length of the file
  flen = as.numeric(strsplit(system(paste("wc -l ",filename,".tmp",sep=""), intern=T),"[[:space:]]")[[1]][1])

  gPos = vector("numeric",flen)
  methProp = vector("numeric",flen)
  index = 1

  ##remove header lines
  a = readLines(con=f,n=1)
  a = readLines(con=f,n=1)

  ##go line by line and adjust the appropriate bins
  for(i in 1:(flen-2)){
    a = readLines(con=f,n=1)
    a = strsplit(a,"\t")[[1]]
    gPos[index] = as.numeric(a[1])
    methProp[index] = as.numeric(a[2])/as.numeric(a[3])
    index = index + 1
  }
  gPos = gPos[1:(index-1)]
  methProp = methProp[1:(index-1)]
  system(paste("rm ",filename,".tmp",sep=""))

  ## convert to IRanges
  meth = IRanges(start=gPos,end=gPos)

  ##figure out which bases fall in which bins
  matches <- as.matrix(findOverlaps(bins,meth))

  ##add the percentage of Cs methylated back in to the number of GC bases
  for(i in 1:length(matches[,1])){
    gcBases[matches[i,1]] = gcBases[matches[i,1]] + methProp[matches[i,2]]
  }
  readBin$gc = gcBases/effectiveLength

  return(readBin)
}



#--------------------------------------------------
## match the GC content values up with the appropriate
## genomic windows

binGC <- function(rdo, chr, readlength){
  binNum = length(rdo@chrs[[chr]][,1])
  len = getChrLength(chr,rdo@params$entrypoints)

  ## read in all gc windows
  gc = scan(gzfile(paste(rdo@params$annotationDirectory,"/gcWinds/",chr,".gc.gz",sep="")), what=0, quiet=TRUE)

  a <- 1:binNum
  numPerBin=rdo@binParams$binSize/rdo@params$gcWindowSize

  getGcMean <- function(x){
    wind = gc[(((x-1)*(numPerBin))+1):(x*(numPerBin))]
    mn = mean(wind, na.rm=TRUE)
    if(is.nan(mn)){
      return(NA)
    } else{
      return(mn)
    }
  }
  return(listify(data.frame(gc=sapply(a,getGcMean)),chr))
}


#--------------------------------------------------
## match the mapability content values up with the appropriate
## genomic windows
binMap <- function(rdo, chr, readlength){
  binNum = length(rdo@chrs[[chr]][,1])
  len = getChrLength(chr,rdo@entrypoints)

  ## read in all mapability windows
  maps = scan(gzfile(paste(rdo@params$annotationDirectory,"/mapability/",chr,".dat.gz",sep="")), what=0, quiet=TRUE)

  a <- 1:binNum
  numPerBin=rdo@binParams$binSize/rdo@params$gcWindowSize

  getMapMean <- function(x){
    wind = maps[(((x-1)*(numPerBin))+1):(x*(numPerBin))]
    mn = mean(wind, na.rm=TRUE)
    if(is.nan(mn)){
      return(NA)
    } else{
      return(mn)
    }
  }
  return(listify(data.frame(map=sapply(a,getMapMean)),chr))
}


##--------------------------------------------------
## remove the specified percentage of outliers from
## the high and low ends to avoid overfitting
##
stripOutliers <- function(nonNAbins,outlierPercentage){
  num=round(length(nonNAbins$avgreads)*outlierPercentage)
  top = sort(nonNAbins$avgreads)[length(nonNAbins$avgreads)-num]
  btm = sort(nonNAbins$avgreads)[1+num]
  rmlist = which(nonNAbins$avgreads >= top)
  rmlist = append(rmlist,which(nonNAbins$avgreads <= btm))
  for(i in 1:length(rmlist)){
    nonNAbins$avgreads[rmlist[i]] = NA
  }

  return(list(nonNAbins,rmlist))
}


##--------------------------------------------------
## restore the outliers, setting their corrected
## value to the mean of the windows around them
returnOutliers <- function(results,rmlist,preOutliers){

  ##first, group adjacent windows that were removed
  counter = 1
  segs = vector("list")
  rmlist = sort(rmlist)
  segs[[1]] = rmlist[1]

  for(i in 2:length(rmlist)){
    ##adjacent windows, merge
    if(rmlist[i]-1 == segs[[counter]][length(segs[[counter]])]){
      segs[[counter]] = append(segs[[counter]],rmlist[i])
    } else { #just add
      counter = counter + 1
      segs[[counter]] = rmlist[i]
    }
  }

  ##now, put them back in
  for(i in 1:length(segs)){
    ##edge case: beginning
    if(1 %in% segs[[i]]){
      ##can't avg, get the  of the closest non-outlier values (high)
      hival = preOutliers[(which(preOutliers$val == preOutliers$val[max(segs[[i]])])+1),]$val
      ##get adjustment for that nearby window
      val = results[which(results$val==hival),]$adj

      for(j in 1:length(segs[[i]])){
        results = rbind(results,data.frame(adj=val,val=preOutliers$val[segs[[i]][j]]))
      }

      ##other edge case: end
    } else if(max(segs[[i]]) > length(preOutliers$reads)){

      ##can't avg, get the content of the closest non-outlier values (low)
      loval = preOutliers[(which(preOutliers$val == preOutliers$val[min(segs[[i]])])-1),]$val
      ##get adjustment for that nearby window
      val = results[which(results$val==loval),]$adj

      for(j in 1:length(segs[[i]])){
        results = rbind(results,data.frame(adj=val,val=preOutliers$val[segs[[i]][j]]))
      }

      ##in the middle, average
    } else{
      ##get the contents of the closest non-outlier values (high and low)
      loval = preOutliers[(which(preOutliers$val == preOutliers$val[min(segs[[i]])])-1),]$val
      hival = preOutliers[(which(preOutliers$val == preOutliers$val[max(segs[[i]])])+1),]$val

      ##now, average the adjustments for those two nearby windows
      val = mean(c(results[which(results$val==hival),]$adj, results[which(results$val==loval),]$adj))

      for(j in 1:length(segs[[i]])){
        results = rbind(results,data.frame(adj=val,val=preOutliers$val[segs[[i]][j]]))
      }
    }
  }
  return(results)
}

##--------------------------------------------------
## correct for bias
##
loessCorrect <- function(rdo, libNum, readLength, libName, outlierPercentage=0.01, corrResolution=0.001, type="gc"){


  binSize=rdo@binParams$binSize
  ##count the number of reads in each percentage bin
  corrBins = c();
  if(type=="gc"){
    chr = NULL
    corrBins = foreach(chr=rdo@entrypoints$chr, .combine="combinePercBins") %dopar% {    
      makeCorrBins(rdo@chrs[[chr]],libNum,corrResolution,chr,type)
    }
  } else { #use median
    corrBins = makeCorrBinsMedian(rdo,libNum,corrResolution,type)
  }

  ##strip out those without any reads
  nonNAbins=na.omit(corrBins)
  ##remove outliers if specified
  rmlist = NULL
  preOutliers=NULL
  if(outlierPercentage > 0){
    preOutliers=nonNAbins
    tmp = stripOutliers(nonNAbins, outlierPercentage)
    nonNAbins = tmp[[1]]
    rmlist = tmp[[2]]
  }
  nonNAbins = na.omit(nonNAbins)
  ##now, do the loess correction
  val = nonNAbins$val
  reads = nonNAbins$avgreads
  numreads = nonNAbins$numreads
  reads.loess <- loess(reads ~ val, span=0.75, data.frame(reads=reads, val=val))
  
  ## do adjustment
  bmed=balancedCenter(reads.loess$fitted,numreads)
#  bmed=rdo@binParams$med
  reads.adj = reads.loess$fitted-bmed
  reads.fix = reads-reads.adj
  reads.resLoess = loess(reads.fix ~ val, span=0.75, data.frame(reads.fix=reads.fix, val=val))
  results = data.frame(adj=reads.adj,val=val)

  ##plot the results of the fit
  ##type=gc|map applies here
  if(type=="map"){
    plotMapLoess(bmed, val, reads, reads.loess$fitted, reads.fix, reads.resLoess$fitted, rdo, libNum, readLength, libName)
  } else if(type=="gc"){
    plotGcLoess(bmed, val, reads, reads.loess$fitted, reads.fix, reads.resLoess$fitted, rdo, libNum, readLength, libName)
  }

  if(outlierPercentage > 0 & length(rmlist) > 0){
    results = returnOutliers(results,rmlist,preOutliers)
  }

  ## restore NA bins we removed at the beginning
  x = 1:length(corrBins$val)
  restoreNAs <- function(z){
    if(!(is.na(corrBins$avgreads[z]))){
      return(results[which(results$val==corrBins$val[z]),]$adj)
    } else {
      return(NA)
    }
  }
  adjustments = sapply(x,restoreNAs)

  ##finally, apply the adjustment to the read bins
  adjBins = foreach(chr=rdo@entrypoints$chr, .combine="combineBins") %dopar% {
    doCorrection(rdo@chrs[[chr]], libNum, corrResolution, adjustments, chr, type)
  }
  closeAllConnections()
  return(adjBins)
}




##---------------------------------------------------------
## choose a median value such that the adjustments
## don't have any effect on the overall median of the
## data set.  This makes sense because we don't account
## for mapability or gc in the inital model, and any adjustments
## to this med will affect every window equally, such that
## the relative copy-number is preserved.
balancedCenter <- function(pos, numReads){
  center = 0
  netChange = sum((pos-center)*numReads)
  while(netChange > 0){
    center = center + 1
    netChange = sum((pos-center)*numReads)
  }

  #now zero in on zero (within 0.1%)
  margin = sum(numReads)*0.001
  adj = 0.5
  center = (center-1)+adj
  netChange = sum((pos-center)*numReads)
  while((abs(netChange) > margin) & (adj > 0)){
    adj = adj/2
    if(netChange > 0){
      center = center + adj
    }else{
      center = center - adj
    }
    netChange = sum((pos-center)*numReads)
  }
  return(center)
}



##---------------------------------------------------------
## apply the calculated correction to the read depths
##
doCorrection <- function(bin, libNum, corrResolution, theAdj, chr, type){

  doAdjust <- function(x){
    #print(x[libNum])
    if( (!(is.na(x[libNum]))) & (!(is.na(x[type])))){
      ##don't correct windows with zero reads
      if(!(x[libNum]==0)){
        x[libNum] <- x[libNum] - theAdj[as.numeric(round(x[type]/corrResolution)+1)]
        ##for the rare case that a bin ends up below zero, set to just
        ## above zero.  (can't be zero, because it has at least one read)
        if(x[libNum] < 0){
          x[libNum] = 0
        }
      }
    }
    return(x)
  }

  ## do the adjustment, then convert the matrix to a dataframe, and return
  ## the adjusted read depth in listified format
  z=data.frame(t(as.matrix(apply(bin,1,doAdjust))))[,libNum]
  return(listify(z,chr))
}





##---------------------------------------------------------
## Combine output of parallel calcs
##
combinePercBins <- function(a,b){

  avgWNAs <- function(x){
    if(is.na(x[1])){
      if(is.na(x[2])){
        return(NA)
      }else{
        return(x[2])
      }
    }else{
      if(is.na(x[2])){
        return(x[1])
      }else{
        return(mean(c(x[1],x[2])))
      }
    }
  }
  df <- data.frame(a$avgreads,b$avgreads)
  a$avgreads = apply(df,1,avgWNAs)

  addWNAs <- function(x){
    if(is.na(x[1])){
      if(is.na(x[2])){
        return(NA)
      }else{
        return(x[2])
      }
    }else{
      if(is.na(x[2])){
        return(x[1])
      }else{
        return(x[1] + x[2])
      }
    }
  }
  df <- data.frame(a$numreads,b$numreads)
  a$numreads = apply(df,1,addWNAs)

  return(a)
}


##---------------------------------------------------------
## calculate the number of reads in each percentage bin
##
makeCorrBins <- function(bin,libNum,windSize,chr,type){
  ##create vectors for percentages and number of reads
  myBin <- c()
  zmean <- c()
  zcnt <- c()
  for(i in 1:((1/windSize)+1)){
    myBin[[i]] <- ((i*windSize)-windSize)
    zmean[[i]] <- NA
    zcnt[[i]] <- 0
  }

  names(bin)[libNum] = "count"
  
  bin = cbind(bin,zbin=(round(bin[[type]]/windSize)+1))
  for(i in 1:length(myBin)){
    zmean[[i]] = mean(bin[bin$zbin == i,]$count, na.rm=T)
    zcnt[[i]] = sum(bin[bin$zbin == i,]$count, na.rm=T)
  }
  return(data.frame(val=myBin, avgreads=zmean, numreads=zcnt))

}

##---------------------------------------------------------
## calculate the number of reads in each percentage bin
##
makeCorrBinsMedian <- function(rdo,libNum,windSize,type){
  bins = makeMapDf(rdo,libNum)
  
  ##create vectors for percentages and number of reads
  myBin <- c()
  zmed <- c()
  zcnt <- c()
  for(i in 1:((1/windSize)+1)){
    myBin[[i]] <- ((i*windSize)-windSize)
    zmed[[i]] <- NA
    zcnt[[i]] <- 0
  }

  bins = cbind(bins,zbin=(round(bins[[type]]/windSize)+1))

  for(i in 1:length(myBin)){
    zmed[[i]] = median(bins[bins$zbin == i,]$count, na.rm=T)
    zcnt[[i]] = sum(bins[bins$zbin == i,]$count, na.rm=T)
  }
  return(data.frame(val=myBin,avgreads=zmed, numreads=zcnt))
}

##----------------------------------------------------------------------------------
## remove counts in that overlap the gaps specified in a file
##
removeGapCounts <- function(rdo,gapsFile){
  gaps = read.table(gapsFile);
  names(gaps) = c("chr","st","sp")

  binSize = rdo@params$binSize

  #this could be vectorized later for speed
  for(i in 1:length(gaps[,1])){
    windst = round(gaps[i,2]/binSize)*binSize
    windsp = (round(gaps[i,3]/binSize)*binSize)+1

    cols = grep("^rd.",names(rdo@chrs[[1]]))
    rdo@chrs[[gaps[i,1]]][windst:windsp,cols]
  }
}


##----------------------------------------------------------------------------------
## remove counts in that overlap the gaps specified in a file
##
coverageFilter <- function(rdo,gapsFile){
  gaps = read.table(gapsFile);
  names(gaps) = c("chr","st","sp")

  binSize = rdo@params$binSize

  #this could be vectorized later for speed
  for(i in 1:length(gaps[,1])){
    windst = round(gaps[i,2]/binSize)*binSize
    windsp = (round(gaps[i,3]/binSize)*binSize)+1

    cols = grep("^rd.",names(rdo@chrs[[1]]))
    rdo@chrs[[gaps[i,1]]][windst:windsp,cols]
  }
}
