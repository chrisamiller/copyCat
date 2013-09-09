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
makeDfLog <-function(binList,params){
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
makeDf <-function(binList,params){
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
## make a dataframe from a paired set of samples with chr, pos, log2 ratio
##
makeDfLogPaired <- function(nrm,tum){
  ##create a merged data frame with windows common to both samples
  dftum = makeDf(tum@chrs,tum@params)
  dfnrm = makeDf(nrm@chrs,nrm@params)
  counts=merge(dftum,dfnrm,by=c("chr","pos"))
  counts = counts[with(counts, order(chr,pos)), ]
  
  counts$score.x = counts$score.x/tum@binParams$med
  counts$score.y = counts$score.y/nrm@binParams$med
  
  df = data.frame(chr=counts$chr,pos=counts$pos,score=log2(counts$score.x/counts$score.y))
}


##---------------------------------------------------
## make a dataframe from a paired set of samples with chr, pos, estimatedCn
##
makeDfPaired <- function(nrm,tum){
  ##create a merged data frame with windows common to both samples
  dftum = makeDf(tum@chrs,tum@params)
  dfnrm = makeDf(nrm@chrs,nrm@params)
  counts=merge(dftum,dfnrm,by=c("chr","pos"))
  counts = counts[with(counts, order(chr,pos)), ]
  
  counts$score.x = counts$score.x/tum@binParams$med
  counts$score.y = counts$score.y/nrm@binParams$med
  
  df = data.frame(chr=counts$chr,pos=counts$pos,score=(counts$score.x/counts$score.y)*2)
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
  #figure out filename
  if(is.null(filename)){
    suffix=".dat";
    if(cnvHmmFormat){
      suffix = ".cnvhmm";
    }
    
    if(!(is.null(rdo@params$prefix))){
      filename=paste(rdo@params$outputDirectory,"/",rdo@params$prefix,".rd.bins.",suffix,sep="")
    } else {
      filename=paste(rdo@params$outputDirectory,"/rd.bins.",suffix,sep="")
    }
  }
   
  #add cnvhmm header if necessary
  appendFile=FALSE
  if(cnvHmmFormat){
    appendFile=TRUE
    ## #WholeGenome_Median:4639
    ## #2xReadCount:4639.00 in 10000 bp
    a1 = paste("#WholeGenome_Median:",rdo@binParams$med,sep="")     
    b1 = paste("#2xReadCount:",rdo@binParams$med," in ",rdo@binParams$binSize," bp",sep="")
    c1 = "CHR\tPOS\tReadCount\tCopyNumber"
    write(c(a1,b1,c1), file=filename)
  } else
    
  df = makeDf(rdo@chrs,rdo@binParams)
  df[,2]=df[,2]-(rdo@params$binSize/2)

  if(cnvHmmFormat){
    cn = (df[,3]/rdo@binParams$med)*2
    df = cbind(df,cn)
  }
  
  ##prevent scientific notation
  options(scipen=999)
  write.table(df, file=filename, row.names=F, col.names=F, quote=F, sep="\t", append=appendFile)
}

writePairedBins <- function(nrm,tum){
  df = makeDfPaired(nrm,tum);
  options(scipen=999)
  filename=paste(tum@params$outputDirectory,"/rd.bins.dat",sep="") 
  write.table(df, file=filename, row.names=F, col.names=F, quote=F, sep="\t")
}

writeThresholds <- function(rdo){
  a1=data.frame(gainThresh=(rdo@binParams$gainThresh/rdo@binParams$med)*2)
  a2=data.frame(lossThresh=(rdo@binParams$lossThresh/rdo@binParams$med)*2)
  a3=data.frame(binSize=rdo@binParams$binSize)
  options(scipen=999)
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


##-----------------------------------------------------------
## dumps the parameters to a file
##
dumpParams <- function(rdo){
  write.table(t(rdo@params),paste(rdo@params$outputDirectory,"/params.",rdo@params$prefix,".txt",sep=""),sep="\t",quote=F,col.names=F)
  write.table(t(rdo@binParams),paste(rdo@params$outputDirectory,"/params.",rdo@params$prefix,".txt",sep=""),sep="\t",quote=F,col.names=F,append=T)
  write.table(rdo@readInfo,paste(rdo@params$outputDirectory,"/libraries.",rdo@params$prefix,".txt",sep=""),sep="\t",quote=F,col.names=T)
}


## ##-----------------------------------------------------------
## ## Author: Kevin Wrightnwith some ideas from Andy Liaw
## ## http://tolstoy.newcastle.edu.au/R/help/04/07/1076.html

## ## x: A data.frame
## ## by: A one-sided formula using + for ascending and - for descending
## ##     Sorting is left to right in the formula

## ## Usage is:
## ## library(nlme);
## ## data(Oats)
## ## sort(Oats, by= ~nitro-Variety)
## sortDataFrame <- function(x, by){
 
##     if(by[[1]] != "~")
##         stop("Argument 'by' must be a one-sided formula.")
 
##     # Make the formula into character and remove spaces
##     formc <- as.character(by[2]) 
##     formc <- gsub(" ", "", formc) 
##     # If the first character is not + or -, add +
##     if(!is.element(substring(formc, 1, 1), c("+", "-")))
##         formc <- paste("+", formc, sep = "")
 
##     # Extract the variables from the formula
##     vars <- unlist(strsplit(formc, "[\\+\\-]"))    
##     vars <- vars[vars != ""] # Remove any extra "" terms
 
##     # Build a list of arguments to pass to "order" function
##     calllist <- list()
##     pos <- 1 # Position of + or -
##     for(i in 1:length(vars)){
##         varsign <- substring(formc, pos, pos)
##         pos <- pos + 1 + nchar(vars[i])
##         if(is.factor(x[, vars[i]])){
##             if(varsign == "-") {
##                 calllist[[i]] <- -rank(x[, vars[i]])
##             } else {
##                 calllist[[i]] <- rank(x[, vars[i]])
##             }
##         } else {
##             if(varsign == "-") {
##                 calllist[[i]] <- -x[, vars[i]]
##             } else {
##                 calllist[[i]] <- x[,vars[i]]
##             }
##         }
##     }
##     return(x[do.call("order", calllist), ])
## }


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

