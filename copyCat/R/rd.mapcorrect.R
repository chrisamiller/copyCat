##--------------------------------------------------
## correct for mapability
##
rd.mapCorrect <- function(rdo, minMapability=0.25){

  if(verbose){
    cat("Correcting for Mapability:  ",date(),"\n")
  }

  ##split the jobs
  mcoptions <- list(preschedule = FALSE)

  bins = foreach(chr=rdo@entrypoints$chr, .combine="combineBins", .options.multicore=mcoptions) %dopar% {
    doMapCorrect(rdo@chrs[[chr]], chr, rdo@binParams, rdo@entrypoints, rdo@params, minMapability)
  }
  
  if(verbose){
    cat("Done:  ",date(),"\n")
  }
  gc()

  ##combine the results
  foreach(chr=rdo@entrypoints$chr) %do% {
    rdo@chrs[[chr]] = bins[[chr]]
  }

  ## if(verbose){
  ##   plotDist(rdo, filename="output/dist.postMapCor.pdf")  
  ## }

  return(rdo)
}



##--------------------------------------------------
## do the mapability correction
##
doMapCorrect <- function(bins, chr, binParams, entrypoints, params, minMapability){
  if(verbose){    
    cat("starting:",chr,"\n")
  }
  ##first, read in the vector of mapability for this chromosome
  len = getChrLength(chr, entrypoints)
  ## read in all mapability windows
  mapVec = scan(gzfile(paste(params$annotationDirectory,"/mapability/",chr,".dat.gz",sep="")), what=0, quiet=TRUE)
  #a <- 1:ceiling(len/binParams$binSize)
  a <- 1:length(bins[,1])

  numPerBin=binParams$binSize/params$gcWindowSize

  getMapMean <- function(x){    
    wind = mapVec[(((x-1)*(numPerBin))+1):(x*(numPerBin))]
    mn = mean(wind, na.rm=TRUE)
    if(is.nan(mn)){
      return(NA)
    } else{
      return(mn)
    }
  }
  bins = cbind(bins,map=sapply(a,getMapMean))
  for(i in 1:params$numLibs){
    bins[,i] = bins[,i] * 1/bins$map
    bins[which(bins$map < minMapability),i] = NA
  }
  return(listify(bins,chr))
}

