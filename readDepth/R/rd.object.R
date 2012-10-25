library('methods')
library('foreach')
library('doMC')
registerDoMC()
library('IRanges')
library('digest')

##-------------------------------------------------
## Set up an object to hold the data, then fill
## it when initialized
##
initRdClass <- function(){
  setClass("rdObject", representation(params="data.frame", binParams="data.frame", entrypoints="data.frame", readInfo="data.frame", chrs="list"))
  setMethod("initialize", "rdObject",
          function(.Object){
            .Object@params=readParameters()
            .Object@entrypoints=readEntrypoints(.Object@params$annotationDirectory)
            .Object@entrypoints=addMapability(.Object@entrypoints,.Object@params$annotationDirectory)
            
            ##only model bins if binSize isn't specified
            if (is.na(.Object@params$binSize)){
              
              if (.Object@params$inputType == "bed"){ 
                verifyFiles(.Object@entrypoints$chr,.Object@params)
              }              
              .Object@readInfo=getReadInfo(.Object@params, .Object@entrypoints)
              .Object@binParams=binParams(.Object@params, .Object@entrypoints, .Object@readInfo)
            }
            
            closeAllConnections()
            return(.Object)              
          })
}

##-------------------------------------------------
## calculate the appropriate bin size,
## gain/loss thresholds, number of chromosomes, etc
##
binParams <-function(params,entrypoints,readInfo){

  ## count the total number of reads
  numReads = getNumReads(readInfo)

  ## get genome size from entrypoints file, adjust by mappability estimates
  genomeSize = sum(entrypoints$length)
  mapPerc=sum(entrypoints$mapPerc*entrypoints$length)/sum(entrypoints$length)
  effectiveGenomeSize = genomeSize * mapPerc

  if(verbose){
    cat(numReads," total reads\n")
    cat("genome size:",genomeSize,"\n")
    cat("genome mapability percentage:",mapPerc,"\n")
    cat("effectiveGenomeSize:",effectiveGenomeSize,"\n")
  }


  ## median value has to be adjusted if we have chromosomes with single ploidy
  ## or expected copy number alterations
  ploidyPerc = ploidyPercentages(effectiveGenomeSize,entrypoints,params)

  if(verbose){
    cat("expect ",
        ploidyPerc$haploidPerc*100,"% haploid,",
        ploidyPerc$diploidPerc*100,"% diploid,",
        ploidyPerc$triploidPerc*100,"triploid\n")
  }

  ## est. coverage of genome by reads
  coverage = numReads * params$readLength / effectiveGenomeSize
  if(verbose){
    cat("approx.",coverage," X coverage of mappable genome \n")
  }

  ## calculate window size based on triploid peak, since it always
  ## produces larger (more conservative) windows
  ploidy = 3
  medAdj = 1
  if(verbose){
    if("medAdjustment" %in% names(params)){
      medAdj = params$medAdjustment
    }
  }

  pTrip <- calcWindParams(numReads=numReads,
                          fdr=params$fdr,
                          genomeSize=effectiveGenomeSize,
                          oDisp=params$overDispersion,
                          ploidy=ploidy,
                          minSize=params$gcWindowSize,
                          ploidyPerc=ploidyPerc,
                          medAdj=medAdj)


  binSize = pTrip$binSize
  binSize = round(binSize/100)*100
  med=pTrip$med
  if(verbose){
    cat("expected mean: ",numReads/(effectiveGenomeSize/binSize),"\n")
    cat("adjusted mean: ",med,"\n")
  }

  ## calculate separation peak for haploid peak
  pHap = fdrRate(binSize, 1, effectiveGenomeSize, numReads, params$overDispersion, ploidyPerc, medAdj)

##   ##plot the output for later review
##   pdf("output/cnSeparation.pdf")
##   plotWindows(pTrip$binSize, effectiveGenomeSize, pHap$div, pTrip$div, params$fdr, numReads, params$overDispersion, ploidyPerc, med)
##   dev.off()

  return(data.frame(binSize=binSize, lossThresh=pHap$div, gainThresh=pTrip$div, med=med, hapPerc=ploidyPerc$haploidPerc, dipPerc=ploidyPerc$diploidPerc, tripPerc=ploidyPerc$triploidPerc))
}



##-------------------------------------------------
## read the parameters from the params file
##
readParameters <- function(){
  ## this gets a little complicated because we want a data
  ## frame that stores both numbers and strings
  if(!(exists("PARAMSFILE"))){
    stop("PARAMSFILE not defined. Please define it before creating the object (PARAMSFILE <<- \"/path/to/paramsfile\")")
  }
  
  params=as.data.frame(t(read.table(PARAMSFILE,row.names=1)))

  ##now convert params to character or numeric, depending
  getParam <- function(df, col){
    return(as.character(df[col][,1]))
  }
  putParam <- function(df, col, val){
    df[col][,1]=val
    return(df)
  }
  for(i in 1:length(params[1,])){
    name = names(params)[i]

    ##these params need to be converted to numbers, the rest stay strings
    if ((name == "readLength") |
        (name == "fdr") |
        (name == "overDispersion") |
        (name == "gcWindowSize") |
        (name == "percCNGain") |
        (name == "percCNLoss") |
        (name == "chunkSize") |
        (name == "maxCores") |
        (name == "readCores") |
        (name == "binSize")){
      val = as.numeric(as.character(getParam(params,name)))
      params = putParam(params,name,val)
    } else {
      val = as.character(getParam(params,name))
      params = putParam(params,name,val)
    }
    
  }

  checkParam <- function(param){
    if(is.na(params[param])){
      stop(paste("ERROR: requred parameter",param,"is missing from params file"))
    }
  }

  ##these params always need to be here
  checkParam("inputType")
  checkParam("readLength")
  checkParam("gcWindowSize")
  checkParam("annotationDirectory")
  checkParam("outputDirectory")

  ##params required for window input
  if(params$inputType == "bins"){
    checkParam("binFile")
    checkParam("binSize")
  ##params required for bam input
  } else if(params$inputType == "bam"){
    checkParam("bamFile")
    if (!(is.na(params$binSize))){
      checkParam("fdr")
      checkParam("overDispersion")
      checkParam("percCNGain")
      checkParam("percCNLoss")
    }
  ##params required for bed input
  } else if(params$inputType == "bed"){
    checkParam("bedDirectory")
    if (!(is.na(params$binSize))){
      checkParam("fdr")
      checkParam("overDispersion")
      checkParam("percCNGain")
      checkParam("percCNLoss")
    }
  }

  if (is.na(params$verbose)){
    verbose <<- FALSE
  } else {
    #there's a better way to do this, but meh.
    if (params$verbose == "1"){
      verbose <<- TRUE
    } else {
      z = tolower(params$verbose)
      if ((z == "true") |
          (z == "t")){
        verbose <<- TRUE           
      }
    }
  }

  ##set multicore options if specified
  if(!(is.na(params$maxCores))){
    options(cores = params$cores)
  }

  ## default value - may change later if we're using bins
  ## as input
  params$numLibs = 1;
  return(params)
}



##--------------------------------------------------
## calculate adjusted genome size, based on the fact that we
## may have haploid chromosomes and/or expected CN alterations
ploidyPercentages <- function(effectiveGenomeSize,ents,params){
  ##first, get the coverage that come from diploid chrs
  diploidPerc = sum((ents$length*ents$mapPerc)[which(ents$ploidy==2)])/effectiveGenomeSize
  diploidPerc = diploidPerc - params$percCNLoss
  diploidPerc = diploidPerc - params$percCNGain

  ##coverage from haploid chrs
  haploidPerc = sum((ents$length*ents$mapPerc)[which(ents$ploidy==1)])/effectiveGenomeSize
  haploidPerc = haploidPerc + params$percCNLoss

  return(data.frame(haploidPerc=haploidPerc,
                    diploidPerc=diploidPerc,
                    triploidPerc=params$percCNGain))
}


##--------------------------------------------------
## returns the percentage of genome covered by
## mappable sequence
##
## if a coverage total file exists, use its value
## else, sum the lengths of the annotations, and
## create the cov total file
##
addMapability <-function(entrypoints,annoDir){
  ##first, we need the mappable regions
  mapTotalFileName = paste(annoDir,"/mapability/totalMappablePerc",sep="")
  mapDir = paste(annoDir,"/mapability/",sep="")
  mapTots = 0
  #default is 100% mapability
  entrypoints = cbind(entrypoints,mapPerc=rep(1,length(entrypoints$chr)))
  tmp=NULL
  #file doesn't exist - create it
  if(!(file.exists(mapTotalFileName))){
    sumMaps <- function(filename){
      a=scan(gzfile(filename),what=0,quiet=TRUE)
      return(sum(a)/length(a))
    }

    tmp=foreach(e=entrypoints$chr, .combine="append") %do%{
      c(e,sumMaps(paste(mapDir,e,".dat.gz",sep="")))
    }
    ##now, place map perc in the appropriate part of the table
    for(i in seq(1,length(tmp),2)){
      entrypoints[which(entrypoints$chr==tmp[i]),]$mapPerc = tmp[i+1]
    }
    write.table(entrypoints[,c("chr","mapPerc")],file=mapTotalFileName,sep="\t",quote=F,row.names=F,col.names=F)
    closeAllConnections()
  }

  tmp = scan(mapTotalFileName,what="",quiet=TRUE)
  for(i in seq(1,length(tmp),2)){
    mapp=as.numeric(tmp[i+1])
    if(tmp[i] %in% entrypoints$chr){
      entrypoints[which(entrypoints$chr==tmp[i]),]$mapPerc = mapp
    }
  }

  closeAllConnections()
  entrypoints$mapPerc= as.numeric(entrypoints$mapPerc)
  return(entrypoints)
}



##--------------------------------------------------
## calculate FDR rate and optimal dividing line
##
fdrRate <- function(windSize,ploidy,genomeSize,numReads,oDisp,ploidyPerc, medAdj){ #nullThis,nullBoth){
  numWinds <- genomeSize/windSize

#  med <- numReads/numWinds
  med <- (numReads/(((genomeSize*ploidyPerc$haploidPerc/2) +
                     (genomeSize*ploidyPerc$triploidPerc*3/2) +
                     (genomeSize*ploidyPerc$diploidPerc)) / windSize))*medAdj

  medAlt <- med*(ploidy/2)

  divFdr=NULL
  if(ploidy < 2){
    divFdr = dividePeaks(medAlt, med, numWinds*ploidyPerc$haploidPerc, numWinds*ploidyPerc$diploidPerc, oDisp)
  } else {
    divFdr = dividePeaks(med, medAlt, numWinds*ploidyPerc$diploidPerc, numWinds*ploidyPerc$triploidPerc, oDisp)
  }
  return(data.frame(fdr=divFdr$fdr, div=divFdr$div, med=med))
}


##--------------------------------------------------
## find the optimal dividing line between the
## two specified peaks. Assumes a poisson
## distribution with overdispersion
##
dividePeaks <- function(amed,bmed,aNum,bNum,oDisp){
  thresholds =(amed+1):(bmed-1)

  mislabeledWinds <- function(thresh){
    low=(1-pnbinom(thresh,size=(amed/oDisp-1),mu=amed))*aNum
    high=(pnbinom(thresh,size=(bmed/oDisp-1),mu=bmed))*bNum
    return(low+high)
  }

  vals=sapply(thresholds,mislabeledWinds)
  lows = which(vals == min(vals))

  ##choose the low point closest
  ##to the halfway point
  diffFromMed = abs((lows+amed)-(amed+bmed/2))
  pos = which(diffFromMed == min(diffFromMed))
  if(length(pos) > 1){
    pos = pos[1]
  }
  return(data.frame(div=lows[pos]+amed,fdr=mislabeledWinds(lows[pos]+amed)/(aNum+bNum)))
}



#-------------------------------------------------
# Calculates the window size that conforms to the
# given FDR rate and the threshold that best
# separates the peaks of ploidy
#
calcWindParams <- function(numReads,fdr,genomeSize,oDisp, ploidy, minSize, ploidyPerc, medAdj, startDiv=100){#nullThis, nullBoth, startDiv=100){
  ## answer has to be this close to the FDR rate (on the lower side)
  tolerance = 0.005
  #starting point for search
  windSize = genomeSize/startDiv
  divider = NULL

  ## first, halve size until we get above FDR threshold
  found <- FALSE
  p <- fdrRate(windSize,ploidy,genomeSize,numReads,oDisp,ploidyPerc, medAdj)

  if(p$fdr > fdr){
    stop("not enough reads to achieve the desired FDR rate")
  }

  while((found == FALSE) & (windSize > minSize)){
    windSize <- round(windSize / 2)
    p <- fdrRate(windSize,ploidy,genomeSize,numReads,oDisp,ploidyPerc, medAdj)

    if(p$fdr > fdr){
      found = TRUE
    }
  }

  if(windSize > minSize){
    ## zero in on a size that's within the desired parameters
    found = FALSE
    adj = windSize/2

    while((found == FALSE) & (windSize > minSize)){
      if(p$fdr > fdr){
        windSize <- round(windSize + adj)
        p <- fdrRate(windSize,ploidy,genomeSize,numReads,oDisp,ploidyPerc, medAdj)

      }else if(p$fdr < (fdr - tolerance)){
        windSize <- round(windSize - adj)
        p <- fdrRate(windSize,ploidy,genomeSize,numReads,oDisp,ploidyPerc, medAdj)

      } else{
        found = TRUE
      }
      adj <- adj/2
    }
  }

  ##if the window size is below the minimum size, have to
  ##recalculate params
  if(windSize < minSize){
    windSize = minSize
  } else {
    ## round to multiple of minSize
    windSize = floor(windSize/minSize)*minSize
  }

  p <- fdrRate(windSize,ploidy,genomeSize,numReads,oDisp,ploidyPerc, medAdj)
  med=p$med
  div=p$div

  return(data.frame(binSize=windSize,div=div,med=med))
}
