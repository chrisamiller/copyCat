##-------------------------------------------------
## gc correction functions
##
rd.gcCorrect <- function(rdo, meth=FALSE, outlierPercentage=0.01){
  if(verbose){
    cat("calculating GC content",date(),"\n")
  }

  #figure out the avg num of reads for each level of GC content
  mcoptions <- list(preschedule = FALSE)
  gcBins = foreach(chr=rdo@entrypoints$chr, .combine="combineBins",.options.multicore=mcoptions) %dopar% {
    binGC(rdo, chr)
  }
  for(i in rdo@entrypoints$chr){
    rdo@chrs[[i]] = cbind(rdo@chrs[[i]],gcBins[[i]])
  }
  closeAllConnections()

  #if bisulfite reads, consider methylation before correcting
  if(meth){
    if(verbose){cat("Correcting GC content for bi-sulfite treatment and methylation:  ",date(),"\n") }    
    foreach(chr=rdo@entrypoints$chr) %dopar%{
      rdo@chrs[[chr]] = methAdjust(rdo,chr)
    } 
  }

  #finally, do the loess correction
  if(verbose){cat("Correcting read depth for GC-content bias:  ",date(),"\n")}

  #for each library
  for(i in 1:rdo@binParams$numLibs){
    if(verbose){cat("correcting library ",i,"\n");}
    name = names(rdo@chrs[[1]][i])
    gcAdj = loessCorrect(rdo,i,outlierPercentage=outlierPercentage,type="gc")
    for(j in rdo@entrypoints$chr){
      rdo@chrs[[j]][[name]]= gcAdj[[j]]
    }
  }

#  if(verbose){
#    plotDist(rdo, filename="output/dist.postGcCor.pdf")  
#    cat("Done",date(),"\n")
#    cat("------------------------\n")
#  }
  return(rdo)
}



##-------------------------------------------------
## mapability correction functions
##
rd.mapCorrect2 <- function(rdo, outlierPercentage=0.01){
  if(verbose){
    cat("calculating mapability content",date(),"\n")
  }

  ##figure out the avg num of reads for each level of mapability
  mcoptions <- list(preschedule = FALSE)
  mapBins = foreach(chr=rdo@entrypoints$chr, .combine="combineBins",.options.multicore=mcoptions) %dopar% {
    binMap(rdo, chr)
  }
  for(i in rdo@entrypoints$chr){
    rdo@chrs[[i]] = cbind(rdo@chrs[[i]],mapBins[[i]])
  }
  closeAllConnections()


  #finally, do the loess correction
  if(verbose){cat("Correcting read depth for mapability bias:  ",date(),"\n")}

  #for each library
  for(i in 1:rdo@binParams$numLibs){
    if(verbose){cat("correcting library ",i,"\n");}
    name = names(rdo@chrs[[1]][i])
    mapAdj = loessCorrect(rdo,i,outlierPercentage=outlierPercentage,type="map")
    for(j in rdo@entrypoints$chr){
      rdo@chrs[[j]][[name]]= mapAdj[[j]]
    }
  }

#  if(verbose){
#    plotDist(rdo, filename="output/dist.postMapCor.pdf")  
#    cat("Done",date(),"\n")
#    cat("------------------------\n")
#  }
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

binGC <- function(rdo, chr){
  binNum = length(rdo@chrs[[chr]][,1])
  len = getChrLength(chr,entrypoints)
  
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
binMap <- function(rdo, chr){
  binNum = length(rdo@chrs[[chr]][,1])
  len = getChrLength(chr,entrypoints)
  
  ## read in all mapability windows
  maps = scan(gzfile(paste(rdo@params$annotationDirectory,"/mappability/",chr,"map.gz",sep="")), what=0, quiet=TRUE)
  
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
## plot the loess correction for gc content
##
plotGcLoess <- function(med,gc,reads,fitted,fixed,fixedfit,params,libNum){
  pdf(paste(outputDirectory,"/loessNorm.lib",libNum,".pdf",sep=""),width=8,height=11)
  par(mfcol=c(2,1))

  ymax <- med*3
  ##plot original
  plot(gc,reads,ylim=c(0,ymax),ylab="mean # reads", main="GC content bias - raw Data", xlab="GC content %")
  points(gc,fitted,col="green",pch=20)
  abline(h=med,col="red")

  ## plot post-adjustment
  plot(gc,fixed,ylim=c(0,ymax),ylab="mean # reads", main="GC content bias - post LOESS correction", xlab="GC content %")
  points(gc,fixedfit,col="green",pch=20)
  abline(h=med,col="red")
  
  dev.off()  
}

##--------------------------------------------------
## plot the loess correction for mapability content
##
plotMapLoess <- function(med,map,reads,fitted,fixed,fixedfit,params,libNum){
  pdf(paste(outputDirectory,"/loessNorm.lib",libNum,".pdf",sep=""),width=8,height=11)
  par(mfcol=c(2,1))

  ymax <- med*3
  ##plot original
  plot(map,reads,ylim=c(0,ymax),ylab="mean # reads", main="Mapability bias - raw Data", xlab="Mapability content %")
  points(map,fitted,col="green",pch=20)
  abline(h=med,col="red")

  ## plot post-adjustment
  plot(map,fixed,ylim=c(0,ymax),ylab="mean # reads", main="Mapability bias - post LOESS correction", xlab="Mapability content %")
  points(map,fixedfit,col="green",pch=20)
  abline(h=med,col="red")
  
  dev.off()  
}


##--------------------------------------------------
## remove the specified percentage of outliers from
## the high and low ends to avoid overfitting
##
stripOutliers <- function(nonNAbins,outlierPercentage){
  num=round(length(nonNAbins$avgreads)*outlierPercentage)
  top = sort(nonNAbins$avgreads)[length(nonNAbins$avgreads)-num]
  btm = sort(nonNAbins$avgreads)[1+num]
  rmlist = which(nonNAbins$avgreads > top)
  rmlist = append(rmlist,which(nonNAbins$avgreads < btm))
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
loessCorrect <- function(rdo, libNum, outlierPercentage=0.01, corrResolution=0.001, type="gc"){

  binSize=rdo@binParams$binSize
  ##count the number of reads in each percentage bin
  corrBins = foreach(chr=rdo@entrypoints$chr, .combine="combinePercBins") %dopar% {    
    makeCorrBins(rdo@chrs[[chr]],libNum,corrResolution,chr)
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
#  if(verbose){
#    plotLoess(bmed, gc, reads, reads.loess$fitted, reads.fix, reads.resLoess$fitted, rdo@params, libNum)
#  }
  
  if(outlierPercentage > 0){
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
    doCorrection(rdo@chrs[[chr]], libNum, corrResolution, adjustments, chr)
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
doCorrection <- function(bin, libNum, corrResolution, theAdj, chr){

  doAdjust <- function(x){    
    #print(x[libNum])
    if( (!(is.na(x[libNum]))) & (!(is.na(x[libNum])))){          
      ##don't correct windows with zero reads
      if(!(x[libNum]==0)){         
        x[libNum] <- x[libNum] - theAdj[round(x["val"]/corrResolution)+1]

        #print(x[libNum])
        ##for the rare case that a bin ends up below zero, set to just
        ## above zero.  (can't be zero, because it has at least one read)
        if(x[libNum] < 0){
          x[libNum] = 1
        }
      }
    }
    return(x)
  }
  
  ## do the adjustment, then convert the matrix to a dataframe, and return
  ## the adjusted read depth in listified format
  z=data.frame(t(as.matrix(apply(bin,1,doAadjust))))[,libNum]
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
makeCorrBins <- function(bin,libNum,windSize,chr){  
  ##create vectors for percentages and number of reads
  myBin <- c()
  numReads <- c()
  numHits <- c()
  for(i in 1:((1/windSize)+1)){
    myBin[[i]] <- ((i*windSize)-windSize)
    numReads[[i]] <- NA
    numHits[[i]] <- 0
  }
  ##fill the vectors
  for(i in 1:length(bin$val)){
    if(!(is.na(bin[,libNum][[i]]))){
      val <- round(bin$val[[i]]/windSize)+1
      if(!is.na(val)){
        ##initialize to zero the first time
        if(is.na(numReads[[val]])){
          numReads[[val]] <- 0
        }
        numReads[[val]] <- numReads[[val]] + bin[,libNum][[i]]
        numHits[[val]] <- numHits[[val]] + 1
      }
    }
  }
  return(data.frame(val=myBin,avgreads=numReads/numHits,numreads=numReads))
}


##--------------------------------------------------
## does what it says on the box
##
replaceNAsWZeros <- function(x){
  if(is.na(x)){
    return(0)
  }
  return(x)
}
