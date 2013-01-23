##-------------------------------------------------
## Common functions used by the copyCat package
##


##-------------------------------------------------
## load libraries

library('methods')
library('foreach')
library('IRanges')
library('doMC')


##-------------------------------------------------
## read in entrypoints file, return a data frame
##
readEntrypoints <- function(annoDir,sex="male"){
  p=read.table(paste(annoDir,"/entrypoints.",sex,sep=""), sep="\t", quote="",
    colClasses=c("character","numeric","numeric"))
  names(p)=c("chr","length","ploidy")  
  return(p)
}


##--------------------------------------------------
## get the number of reads in the bam file from samtools flagstat
##
getNumReads <- function(inputFile,outputDirectory){
  ##first check to see if a flagstat exists next to the bam
  flagstat = paste(inputFile,".flagstat",sep="")
  ##if not, run flagstat
  if(!(file.exists(flagstat))){
    print("running samtools flagstat to get read counts")
    system(paste("samtools flagstat ",inputFile," >",outputDirectory,"/flagstat",sep=""));
    flagstat = paste(outputDirectory,"/flagstat",sep="")
  }
  return(strsplit(readLines(flagstat)[5], " ")[[1]][1])
}



##--------------------------------------------------
##  make sure each entrypoint has a corresponding bed file
##
verifyFiles <- function(chrs,params){
  for(i in 1:length(chrs)){
    filename = paste(params$bedDirectory,"/",chrs[i],".bed",sep="")
    if(!(file.exists(filename))){
      stop("file: '",filename,"' does not exist\n")
    }
  }
}


##-------------------------------------------------
## sum up the already-stored # of reads
##
getNumReads <- function(flist){
  len = sum(as.numeric(flist$len))
  return(len)
}

##--------------------------------------------------
## generate a distribution by modelling the overdispersed
## poisson with the negative binomial distribution
## d = var/mean
##
rpois.od<-function (n, lambda,d=1) {
  if (d==1)
    rpois(n, lambda)
  else
    rnbinom(n, size=(lambda/(d-1)), mu=lambda)
}


##----------------------------------------------------
## return the length or ploidy for a given chr identifier
## 
getChrLength <- function(chr,entries){
  return(entries[which(entries$chr==chr),]$length)
}

getChrPloidy <- function(chr,entries){
  return(entries[which(entries$chr==chr),]$ploidy)
}


##----------------------------------------------------
## create a list of the appropriate size and place
## input in the appropriate bin 
##
listify <- function(input,chr){
  lst = vector("list")
  lst[[chr]] = input
  return(lst)
}


##----------------------------------------------------
## function to combine lists returned from parallel
## intersection procedure
##
combineBins <- function(a,b){ 
  for(i in 1:length(a)){
    name=names(a)[i]
    if(is.null(b$name)){
      b[[name]] = a[[i]]
    }
  }
  return(b)
}
 

##---------------------------------------------------
## convert chromosome number to chr name
##
chrName <- function(num){
  if(num == 23){ num <- "X"}
  if(num == 24){ num <- "Y"}
  
  return(paste("chr",num,sep=""))  
}


##---------------------------------------------------
## sum the total lengths of the annotations in a 
## bed file. This is equal to coverage if annotations
## are non-overlapping
##
bedAnnotationLength <- function(e){
  if(file.exists(e)){
    a=scan(gzfile(e),what=0,quiet=TRUE)
    closeAllConnections()
    return(sum((a[seq(2,(length(a)),2)]-a[seq(1,(length(a)-1),2)]+1)))
  }
  return(0)
}


##--------------------------------------------------
## convert read depth to log2 value,
## based on median
##
logScore <- function(val,med){
  if(is.na(val)){
    return(NA)

  }else if(val<=0){
    return(floor(log2(1/med)))

  }else{
    return(log2(val/med))
  }    
}


##--------------------------------------------------
## get the median read count
##
getMedianReadCount <- function(rdo){
  tmp = c()
  for(chr in rdo@entrypoints$chr){    
    tmp = c(tmp,rdo@chrs[[chr]]$rd)
    rdo@binParams$med
  }
}

##--------------------------------------------------
## merge libraries into a single set of reads
##
mergeLibraries <- function(rdo){
  counts = c();
  ##for each chromosome  
  for(chr in rdo@entrypoints$chr){    
    libNames = names(rdo@chrs[[1]])[grep("^rd.",names(rdo@chrs[[1]]))]
    tmp = rdo@chrs[[chr]][[libNames[1]]];

    if(rdo@binParams$numLibs > 1){
      for(i in 2:rdo@binParams$numLibs){
        name=libNames[i]
        tmp = tmp + rdo@chrs[[chr]][[name]];
      }
    }
    rdo@chrs[[chr]][["rd"]] = tmp;
    counts = c(counts,rdo@chrs[[chr]][["rd"]])
  }

  rdo@binParams$med = median(counts, na.rm=T)
  return(rdo);
}

##---------------------------------------------------
## make a dataframe from our rdo object with chr, pos, log score
##
makedfLog <-function(binList,params){
  #drop the last window from each chr, because it's truncated
  for(i in names(binList)){
    binList[[i]] = binList[[i]][1:(length(binList[[i]][,1])-1), ,drop=FALSE]
  }

  ##convert rd to log2 score, make a big list
  logs <- foreach(i=names(binList), .combine="append") %do% {
    sapply(binList[[i]]$rd,logScore,params$med)
  }

  ##generate positions, put in a big list
  pos <- foreach(i=names(binList), .combine="append") %do% {
    #marker is at beginning of bin,  move to center of bin
    ((0:(length(binList[[i]]$rd)-1))*params$binSize)+(params$binSize/2)
  }

  #generate list of chromomsome markers
  chrom <- foreach(i=names(binList), .combine="append") %do% {
    rep(i,length(binList[[i]]$rd))
  }
  #combine them all into properly formatted data frame
  gd = data.frame(chr=chrom,pos=pos,score=logs)

  return(subset(gd,!is.na(score)))
}

##---------------------------------------------------
## make a dataframe from our rdo object with chr, pos, score
##
makedf <-function(binList,params){
  #drop the last window from each chr, because it's truncated
  for(i in names(binList)){
    binList[[i]] = binList[[i]][1:(length(binList[[i]][,1])-1), ,drop=FALSE]
  }

  ##get scores 
  logs <- foreach(i=names(binList), .combine="append") %do% {
    binList[[i]]$rd
  }

  ##generate positions, put in a big list
  pos <- foreach(i=names(binList), .combine="append") %do% {
    #marker is at beginning of bin,  move to center of bin
    ((0:(length(binList[[i]]$rd)-1))*params$binSize)+(params$binSize/2)
  }

  #generate list of chromomsome markers
  chrom <- foreach(i=names(binList), .combine="append") %do% {
    rep(i,length(binList[[i]]$rd))
  }
  #combine them all into properly formatted data frame
  gd = data.frame(chr=chrom,pos=pos,score=logs)

  return(subset(gd,!is.na(score)))
}


##---------------------------------------------------
## make a dataframe from a particular library with count
## and mapablity
##
makeMapDf <-function(rdo,libNum){
  binList = rdo@chrs
  ##drop the last window from each chr, because it's truncated  
  for(i in names(rdo@chrs)){
    binList[[i]] = rdo@chrs[[i]][1:(length(rdo@chrs[[i]][,1])-1), ,drop=FALSE]
  }

  ##get scores 
  counts <- foreach(i=names(binList), .combine="append") %do% {
    binList[[i]][[libNum]]
  }

  ##get map
  maps <- foreach(i=names(binList), .combine="append") %do% {
    binList[[i]]$map
  }
  
  return(data.frame(count=counts,map=maps))
}

##-----------------------------------------------------------
## some simple output functions
##
writeBins <- function(rdo, filename=NULL, cnvHmmFormat=FALSE){#, includeMapability=FALSE, includeGcContent=FALSE){

  #add cnvhmm header if necessary
  appendFile=FALSE
  if(cnvHmmFormat){
    appendFile=TRUE
    if(is.null(filename)){
      filename=paste(rdo@params$outputDirectory,"/rd.bins.cnvhmm",sep="")
    }
    
    ## #WholeGenome_Median:4639
    ## #2xReadCount:4639.00 in 10000 bp
    a1 = paste("#WholeGenome_Median:",rdo@binParams$med,sep="")     
    b1 = paste("#2xReadCount:",rdo@binParams$med," in ",rdo@binParams$binSize," bp",sep="")
    c1 = "CHR\tPOS\tReadCount\tCopyNumber"
    write(c(a1,b1,c1), file=filename)
  }
    
  if(is.null(filename)){
    filename=paste(rdo@params$outputDirectory,"/rd.bins.dat",sep="")
  }
  
  df = makedf(rdo@chrs,rdo@binParams)
  df[,2]=df[,2]-(rdo@params$binSize/2)

  if(cnvHmmFormat){
    cn = (df[,3]/rdo@binParams$med)*2
    df = cbind(df,cn)
  }
  
  ##prevent scientific notation
  options(scipen=999)
  write.table(df, file=filename, row.names=F, col.names=F, quote=F, sep="\t", append=appendFile)
}


writeThresholds <- function(rdo){
  a1=data.frame(gainThresh=(rdo@binParams$gainThresh/rdo@binParams$med)*2)
  a2=data.frame(lossThresh=(rdo@binParams$lossThresh/rdo@binParams$med)*2)
  a3=data.frame(binSize=rdo@binParams$binSize)
  write.table(t(cbind(a1,a2,a3)),file=paste(rdo@params$outputDirectory,"/thresholds.dat",sep=""),sep="\t",quote=F,row.names=T,col.names=F)
}

writeSegs <- function(segs,rdo,filename="segs.dat"){
  options(scipen=999)
  write.table(segs, file=paste(rdo@params$outputDirectory,"/",filename,sep=""), sep="\t", quote=F, row.names=F, col.names=F)
}

writeAlts <- function(segs,rdo,filename="alts.dat"){
  options(scipen=999)
  write.table(getAlts(segs,rdo), file=paste(rdo@params$outputDirectory,"/",filename,sep=""), sep="\t", quote=F, row.names=F, col.names=F)
}



##------------------------------------------------------------
## Estimate overdispersion in a set of reads
## 
estimateOd <- function(rdo, maxOd=30, histBreaks=10){
  #get the bins together
  bins=2^(makedf(rdo@chrs,rdo@binParams)$score)*rdo@binParams$med

  fracNormal=rdo@binParams$dipPerc
  
  ## loop through different overdispersion values, see which one
  ## fits the data best
  errs = c();
  for(i in 1:maxOd){
    model=rpois.od(length(bins)*fracNormal,median(bins,na.rm=T),i)

    maxBin = max(max(bins),max(model))
    
    real.hist = hist(bins, breaks=seq(0,maxBin+histBreaks,histBreaks), 
      xlim=c(0,median(bins,na.rm=T)*2), col=rgb(1, 0, 0,0.5),
      main=paste("Read depth distribution - Model overdispersion:",i))

    model.hist = hist(model, breaks=seq(0,maxBin+histBreaks,histBreaks),
      xlim=c(0,median(bins,na.rm=T)*2), col=rgb(0, 0, 1,0.5),
      add=T)

    mtext("real=red/model=blue")
    
    # here's the fit test - just a simple difference of counts in each bin
    error = sum(abs(real.hist$counts-model.hist$counts))
    errs = c(errs,error)    
  }
  
  bestod = which(errs == min(errs))

#  pdf(paste(rdo@params$outputDirectory,"/odfit.pdf",sep=""));

  model=rpois.od(length(bins)*fracNormal,median(bins,na.rm=T),bestod)

  real.hist = hist(bins, breaks=seq(0,maxBin+histBreaks,histBreaks), 
    xlim=c(0,median(bins,na.rm=T)*2), col=rgb(1, 0, 0,0.5),
    main=paste("Read depth distribution - Model overdispersion:", bestod))
  
  model.hist = hist(model, breaks=seq(0,maxBin+histBreaks,histBreaks),
    xlim=c(0,median(bins,na.rm=T)*2), col=rgb(0, 0, 1,0.5),
    add=T)

  mtext("real=red/model=blue")

#  dev.off()
  if(bestod >= maxOd){
    return(paste("greater than ",maxOd,sep=""))
  }
  return(bestod)
}




##------------------------------------------------------------
## Code for pulling down the annotation files from google code
## and sticking them in the right place

##TODO - fix this to reflect new directory structure

getAnnotations <- function(readLength, sex, genome="hg18", bs=FALSE, annoDir){
  if(sex!="male" & sex!="female" & sex!="autosomes"){
    print("'sex' must be either \"male\", \"female\", or \"autosomes\"")
    return(0)
  }

  oldloc = getwd()
  dir.create(annoDir, showWarnings=FALSE)
#  setwd(annoDir)
  
  dlAnnFile <- function(url,outfile){
    download.file(url, outfile, quiet=TRUE)
    if(file.exists(outfile)){
      print(paste(url,"downloaded successfully"))
      return(1)
    } else {
      print(paste("Oops! Couldn't fetch ",url))
      print("Either you specified a read length or genome build that we don't have annotation")
      print("files prepared for, or your network connection isn't working.")
      print("Check the documentation at http://code.google.com/p/readdepth/ for more information")
      return(0)
    }
  }


  entryurl=paste("http://readdepth.googlecode.com/files/entrypoints.",genome,".",sex,sep="")
  entryfile=paste(annoDir,"/entrypoints",sep="")    
  
  if(bs){ #bisulfite
    mapurl = paste("http://readdepth.googlecode.com/files/mapability.bs.readLength",readLength,".",genome,".tar",sep="")
    mapfile = paste(annoDir,"/mapability.bs.readLength",readLength,".",genome,".tar",sep="")
    gcurl = paste("http://readdepth.googlecode.com/files/gcWinds.bs.readLength",readLength,".",genome,".tar",sep="")
    gcfile = paste(annoDir,"/gcWinds.bs.readLength",readLength,".",genome,".tar",sep="")    
  } else { #normal
    mapurl = paste("http://readdepth.googlecode.com/files/mapability.readLength",readLength,".",genome,".tar",sep="")
    mapfile = paste(annoDir,"/mapability.readLength",readLength,".",genome,".tar",sep="")
    gcurl = paste("http://readdepth.googlecode.com/files/gcWinds.readLength",readLength,".",genome,".tar",sep="")
    gcfile = paste(annoDir,"/gcWinds.readLength",readLength,".",genome,".tar",sep="")    
  }

  #download the files, untar them
  if(dlAnnFile(mapurl, mapfile)){
    system(paste("tar -xf ",mapfile," -C ",annoDir,"/",sep=""))
  }
  
  if(dlAnnFile(gcurl, gcfile)){
    system(paste("tar -xf ",gcfile," -C ",annoDir,"/",sep=""))
  }

  dlAnnFile(entryurl, entryfile)

  setwd(oldloc)
  ##    system(paste("tar -xf",entryfile,"-C annotations/"))
#  }
}


##-----------------------------------------------------------
chooseAnnotationReadLength <- function(readLength, availableReadLengths){
  diff=sqrt((availableReadLengths-readLength)^2)
  return(availableReadLengths[tail(which(diff == min(diff)),1)])
}


##-----------------------------------------------------------
## takes two rd objects and sums the values
## in their rd columns
addObjectBins <- function(rdo,rdo2){
  counts = c();
  ##for each chromosome  
  for(chr in rdo@entrypoints$chr){    
    rdo@chrs[[chr]][["rd"]] = rdo@chrs[[chr]][["rd"]] + rdo2@chrs[[chr]][["rd"]];
  }
  return(rdo);
}



##-----------------------------------------------------------
## takes two rd object representing a tumor/normal pair and
## creates a file suitable for input to CNVHMM
## This mimics the Bam2CNA functionality
##
writeCnahmmInput <- function(nrm,tum){
  dftum = makedf(tum@chrs,tum@params)
  dfnrm = makedf(nrm@chrs,nrm@params)

  #get rid of the old file if one exists
  if(file.exists(paste(tum@params$outputDirectory,"/tumor.normal.cn",sep=""))){
    unlink(paste(tum@params$outputDirectory,"/tumor.normal.cn",sep=""))
  }

  for(chr in unique(dftum[,1])){
#    tumcnts = dftum[which(dftum[,1] == chr),]
#    nrmcnts = dfnrm[which(dfnrm[,1] == chr),]
#    counts=merge(tumcnts,nrmcnts,by="pos")

#    x = paste("#Chr",chr," running average read count\ttumor:",mean(tumcnts),"\tnormal:",mean(nrmcnts),sep="")
    x = paste("#Chr",chr," running average read count\ttumor:",tum@binParams$med,"\tnormal:",nrm@binParams$med,sep="")
    write.table(x, file=paste(tum@params$outputDirectory,"/tumor.normal.cn",sep=""), append=TRUE,
                quote=FALSE, sep="\t", col.names=FALSE, row.names=FALSE)
  }

  x = "CHR\tPOS\tTUMOR\tNORMAL\tDIFF"
  write.table(x, file=paste(tum@params$outputDirectory,"/tumor.normal.cn",sep=""), append=TRUE,
              quote=FALSE, sep="\t", col.names=FALSE, row.names=FALSE)


  counts=merge(dftum,dfnrm,by=c("chr","pos"))
  counts = sort(counts,by = ~ +chr +pos)  
  ratio = sum(counts$score.x)/sum(counts$score.y)
  
  for(chr in unique(counts$chr)){
    ## when calculating diff, downsample the one with the greater number
    ## of reads so they have comparable coverage
    tumcnts = counts[which(counts$chr == chr),]$score.x
    nrmcnts = counts[which(counts$chr == chr),]$score.y

    
    diffScore = c()
    if(ratio > 1){
      diffScore = ((tumcnts*(1-ratio) - nrmcnts)*2)/mean(nrmcnts)
    } else {
      diffScore = ((tumcnts - nrmcnts*(ratio))*2)/mean(nrmcnts)
    }

    ## a = data.frame(chr=dftum[which(dftum[,1]==chr),1],
    ##            pos=dftum[which(dftum[,1]==chr),2],
    ##            tum=dftum[which(dftum[,1]==chr),3],
    ##            nrm=dfnrm[which(dfnrm[,1]==chr),3])
      
    a = cbind(counts[which(counts$chr == chr),],diff=diffScore)
    write.table(a, file=paste(tum@params$outputDirectory,"/tumor.normal.cn",sep=""), append=TRUE,
                quote=FALSE, sep="\t", col.names=FALSE, row.names=FALSE)
    
  }
}

##-----------------------------------------------------------
## takes one rd object and
## creates a file suitable for input to CNVHMM
## This mimics the Bam2CN functionality
##
writeCnvhmmInput <- function(tum){
  dftum = makedf(tum@chrs,tum@params)

  #get rid of the old file if one exists
  if(file.exists(paste(tum@params$outputDirectory,"/sample.cn",sep=""))){
    unlink(paste(tum@params$outputDirectory,"/sample.cn",sep=""))
  }

  x = paste("#WholeGenome_Median:",tum@binParams$med,sep="")
  write.table(x, file=paste(tum@params$outputDirectory,"/sample.cn",sep=""), append=TRUE,
              quote=FALSE, sep="\t", col.names=FALSE, row.names=FALSE)

#  x = paste("#2xReadCount:",tum@binParams$med," in ",????sep="")
#  write.table(x, file=paste(tum@params$outputDirectory,"/sample.cn",sep=""), append=TRUE,
#              quote=FALSE, sep="\t", col.names=FALSE, row.names=FALSE) 
 
  x = "CHR\tPOS\tReadCount\tCopyNumber"
  write.table(x, file=paste(tum@params$outputDirectory,"/sample.cn",sep=""), append=TRUE,
              quote=FALSE, sep="\t", col.names=FALSE, row.names=FALSE)


  dftum = cbind(dftum$chr, dftum$pos, dftum$score, (dftum$score/tum@binParams$med)*2)
  write.table(dftum, file=paste(tum@params$outputDirectory,"/sample.cn",sep=""), append=TRUE,
                quote=FALSE, sep="\t", col.names=FALSE, row.names=FALSE)    
}


##-----------------------------------------------------------
## dumps the parameters to a file
##
dumpParams <- function(rdo){
  write.table(t(rdo@params),paste(rdo@params$outputDirectory,"/params.",rdo@params$prefix,".txt",sep=""),sep="\t",quote=F,col.names=F)
  write.table(t(rdo@binParams),paste(rdo@params$outputDirectory,"/params.",rdo@params$prefix,".txt",sep=""),sep="\t",quote=F,col.names=F,append=T)
  write.table(rdo@readInfo,paste(rdo@params$outputDirectory,"/libraries.",rdo@params$prefix,".txt",sep=""),sep="\t",quote=F,col.names=T)
}


##-----------------------------------------------------------
## Author: Kevin Wright
## with some ideas from Andy Liaw
## http://tolstoy.newcastle.edu.au/R/help/04/07/1076.html

## x: A data.frame
## by: A one-sided formula using + for ascending and - for descending
##     Sorting is left to right in the formula

## Useage is:
## library(nlme);
## data(Oats)
## sort(Oats, by= ~nitro-Variety)

sort.data.frame <- function(x, by){
 
    if(by[[1]] != "~")
        stop("Argument 'by' must be a one-sided formula.")
 
    # Make the formula into character and remove spaces
    formc <- as.character(by[2]) 
    formc <- gsub(" ", "", formc) 
    # If the first character is not + or -, add +
    if(!is.element(substring(formc, 1, 1), c("+", "-")))
        formc <- paste("+", formc, sep = "")
 
    # Extract the variables from the formula
    vars <- unlist(strsplit(formc, "[\\+\\-]"))    
    vars <- vars[vars != ""] # Remove any extra "" terms
 
    # Build a list of arguments to pass to "order" function
    calllist <- list()
    pos <- 1 # Position of + or -
    for(i in 1:length(vars)){
        varsign <- substring(formc, pos, pos)
        pos <- pos + 1 + nchar(vars[i])
        if(is.factor(x[, vars[i]])){
            if(varsign == "-") {
                calllist[[i]] <- -rank(x[, vars[i]])
            } else {
                calllist[[i]] <- rank(x[, vars[i]])
            }
        } else {
            if(varsign == "-") {
                calllist[[i]] <- -x[, vars[i]]
            } else {
                calllist[[i]] <- x[,vars[i]]
            }
        }
    }
    return(x[do.call("order", calllist), ])
}

##-----------------------------------------------------------
## pull out the subset of reads with a specific read length
##
subsetByReadLength <- function(rdo,length){  
  ##get the names of columns with this read length
  libs = rdo@readInfo[which(rdo@readInfo$readlength==length),]$lib
  tmp=c();for(l in libs){tmp=c(tmp,paste("rd.",l,".",length,sep=""))}
  libs = tmp;
  cols=names(rdo@chrs[[1]]) %in% libs
  #do the subset
  for(i in names(rdo@chrs)){    
    rdo@chrs[[i]] = rdo@chrs[[i]][cols]
  }  
  return(rdo)
}



##-----------------------------------------------------------
## restore reads that have been corrected (in rdo2) back into rdo
##
replaceReadCounts <- function(rdo,rdo2){  
  ##get the names of columns in rdo2
  libs=names(rdo2@chrs[[1]])
  for(i in 1:length(libs)){
    for(chr in names(rdo@chrs)){    
      rdo@chrs[[chr]][which(names(rdo@chrs[[1]]) %in% libs[i])] = rdo2@chrs[[chr]][i]
    }
  }
  return(rdo)
}

##-----------------------------------------------------------
## find the appropriate annotations to use with this data
## if it doesn't exist, warn the user
##
getAnnoDir <- function(annodir, readlength, tolerance=5){
  idealDir = paste(annodir,"/readlength.",readlength,sep="")
  if(file.exists(idealDir)){
    return(idealDir)
  }
  ##else, see what other directories are available
  getlens <- function(dir){
    a = strsplit(dir,".",fixed=T)[[1]]
    return(as.numeric(a[length(a)]))
  }

  dirs = sapply(Sys.glob(paste(annodir,"/readlength*",sep="")), getlens)
  diff = sort(abs(dirs-readlength))
  if(diff[1] > tolerance){    
    print(paste("ERROR: no annotation files exists that match a read length of ",readlength," (+/-",tolerance,") ",sep=""))
    print("consult the documentation for instructions on downloading or creating these files.")
    stop()
  }
  if(verbose){
    print(paste("WARNING: annotations for read length of ",readlength," don't exist",sep=""))
    print(paste("using annotations for read length of ",getlens(names(diff[1]))," which are close enough",sep=""))
  }
  return(names(diff[1]))
}

