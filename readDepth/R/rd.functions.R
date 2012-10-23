##-------------------------------------------------
## Common functions used by the readDepth package
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
readEntrypoints <- function(annoDir){
  p=read.table(paste(annoDir,"/entrypoints",sep=""), sep="\t", quote="",
    colClasses=c("character","numeric","numeric"))
  names(p)=c("chr","length","ploidy")  
  return(p)
}


##--------------------------------------------------
##  Do a wc on each read file to find out how many reads there 
##  are. This is way faster than counting with R. Eventually,
##  we can write a little C function to do this so we don't  
##  require a shell with wc installed (allow windows port)
##
getReadInfo <- function(params, entrypoints){
  
  ## check to see if they've already been counted
  lenFileName=paste(params$readDirectory,"/readInfo",sep="")
  len = 0
  
  if(file.exists(lenFileName)){
    flist = read.table(lenFileName)
    names(flist) = c("len","file")
  } else {
    
    ##count up rads in each bed file
    if(params$inputType == "bed"){
      lenList = foreach(e=Sys.glob(paste(params$bedDirectory,"/*.bed"),sep=""), .combine="append") %do%{
        strsplit(system(paste('wc -l ', e, sep=""), intern=T),"[[:space:]]")[[1]]
      }
      flist = data.frame(len=as.numeric(lenList[seq(1,length(lenList)-1,2)]), file=lenList[seq(2,length(lenList),2)])
      
    ##count reads from single bam, keep those that match our entrypoints
    } else if(params$inputType == "bam"){        
      a = read.table(pipe(paste("samtools view -c -q 1 ", params$bamFile," | cut -f 3 | sort | uniq -c")))
      valid = which(a %in% entrypoints$chr)
      flist = data.frame(len=a[valid,1],file=a[valid,V2])
    }
    write.table(flist,file=lenFileName,sep="\t",quote=FALSE, row.names=FALSE, col.names=FALSE)

  }  
  return(flist)
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

    if(rdo@params$numLibs > 1){
      for(i in 2:rdo@params$numLibs){
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


##----------------------------------------------
## given a set of parameters, plot the peaks
## and thresholds
##
plotWindows <- function(windowSize, genomeSize, divGain, divLoss, fdr, numReads, oDisp, ploidyPerc, med){
  numWinds <- genomeSize/windowSize
  
  ##generate distributions
  d2 <- rpois.od(numWinds*ploidyPerc$diploidPerc,med,oDisp)
  d3 <- rpois.od(numWinds*ploidyPerc$triploidPerc,(med*1.5),oDisp)
  d1 <- rpois.od(numWinds*ploidyPerc$haploidPerc,(med*0.5),oDisp)

  hrange=10*sqrt(med*oDisp)
  hist(d2,breaks=seq(-100000,100000,20),xlim=c(med-hrange,med+hrange),main=paste("Window Size: ",windowSize,sep=""))
  mtext(paste("FDR: ",round(fdr,5),sep=""))
  hist(d3,breaks=seq(-100000,100000,20),add=T,col="red")
  hist(d1,breaks=seq(-100000,100000,20),add=T,col="red")
  abline(v=med,col="blue")
  abline(v=med*(0.5),col="blue")
  abline(v=med*(1.5),col="blue")
  abline(v=divGain,col="green")
  abline(v=divLoss,col="green")
}

##----------------------------------------------
## plot the segments for a given chromosome
##
plotSegs <- function(rdo,segs,chr){

  st = 1
  sp = rdo@entrypoints[which(rdo@entrypoints$chr == chr),]$length
  
  winds = rdo@chrs[[chr]]$rd
  #print(winds)
  binSize = rdo@binParams$binSize  
  pos = seq(binSize/2,(((length(winds)-1)*binSize)-binSize/2),binSize)
  pos = append(pos,sp)
  
  par(mar=c(5, 4, 4, 4) + 0.1)
  plot(pos,winds,ylab="number of reads", xlab="position (bp)")

  abline(h=rdo@binParams$med,col="blue")
  abline(h=rdo@binParams$gainThresh,col="green")
  abline(h=rdo@binParams$lossThresh,col="green")

  asegs = segs[which(segs$chrom == chr),]
  for(i in 1:length(asegs[,1])){
    lines(c(asegs[i,2],asegs[i,3]), c(asegs[i,5],asegs[i,5]), col="red",lwd=3)
  }
  
  par(new=T)
  plot(-10000,-10000,ylim=c(0,(max(winds,na.rm=TRUE)/rdo@binParams$med)), xlim=c(1,sp),axes=F,xlab="", ylab="")
  axis(4, ylim=c(0,max(winds,na.rm=TRUE)/rdo@binParams$med), col="red",col.axis="red")
  mtext("Copy Number",side=4,col="red",line=2.5)

}


##--------------------------------------------------
## plot overlapping histograms
##
plotOverlappingHist <- function(a, b, colors=c("white","gray20","gray50"),
                                breaks=NULL, xlim=NULL, ylim=NULL){
  
  ahist=NULL
  bhist=NULL

  if(!(is.null(breaks))){
    ahist=hist(a,breaks=breaks,plot=F)
    bhist=hist(b,breaks=breaks,plot=F)
  } else {
    ahist=hist(a,plot=F)
    bhist=hist(b,plot=F)

    dist = ahist$breaks[2]-ahist$breaks[1]
    breaks = seq(min(ahist$breaks,bhist$breaks),max(ahist$breaks,bhist$breaks),dist)

    ahist=hist(a,breaks=breaks,plot=F)
    bhist=hist(b,breaks=breaks,plot=F)
  }

  if(is.null(xlim)){
    xlim = c(min(ahist$breaks,bhist$breaks),max(ahist$breaks,bhist$breaks))
  }

  if(is.null(ylim)){
    ylim = c(0,max(ahist$counts,bhist$counts))
  }

  overlap = ahist
  for(i in 1:length(overlap$counts)){
    if(ahist$counts[i] > 0 & bhist$counts[i] > 0){
      overlap$counts[i] = min(ahist$counts[i],bhist$counts[i])
    } else {
      overlap$counts[i] = 0
    }
  }

  plot(ahist, xlim=xlim, ylim=ylim, col=colors[1])
  plot(bhist, xlim=xlim, ylim=ylim, col=colors[2], add=T)
  plot(overlap, xlim=xlim, ylim=ylim, col=colors[3], add=T)
}


##--------------------------------------------------
## plot actual and expected histograms
##
plotDist <- function(rdo,xmax=NULL,filename="output/hist.pdf",windSize=NULL){

  bins=c()
  for(i in 1:length(names(rdo@chrs))){
    bins = append(bins,rdo@chrs[[names(rdo@chrs)[i]]]$rd)
  }

  if(is.null(xmax)){
    xmax=rdo@binParams$gainThresh*2 ##max(bins,na.rm=TRUE) ##sort(bins[round(length(bins)*0.999)])
  }
  
  if(is.null(windSize)){
    windSize = round(xmax/100)
  }

  len=length(bins)
  e=rdo@entrypoints
  breaks=seq(0,xmax+windSize,windSize)
  oDisp=rdo@params$overDispersion

  med=rdo@binParams$med 

  p1 = rpois.od(len*rdo@binParams$hapPerc,med/2,oDisp)
  p2=rpois.od(len*rdo@binParams$dipPerc,med,oDisp)
  p3 = rpois.od(len*rdo@binParams$tripPerc,(med*(3/2)),oDisp)
  
  p=append(append(p1,p2),p3)

  bins = bins[which(bins<xmax & bins>0)]
  p = p[which(p<xmax & p>0)]

  pdf(file=filename)
  plotOverlappingHist(a=bins,b=p,breaks=breaks,xlim=c(0,xmax))
  abline(v=rdo@binParams$gainThresh,col="red")
  abline(v=rdo@binParams$lossThresh,col="red")
  dev.off()
}


##--------------------------------------------------
## plot three actual and expected histograms 
## usually used with raw read depth, post mapability
## correction, and post gc content correction
##
tripDist <- function(xmax=NULL,filename="output/hist.pdf",windSize=20, rdo2, rdo3, rdo4){
  pdf(file=filename, width=12,height=4)
  par(mfcol=c(1,3))
  
  for(rdo in c(rdo2,rdo3,rdo4)){
    bins=c()
    for(i in 1:length(names(rdo@chrs))){
      bins = append(bins,rdo@chrs[[names(rdo@chrs)[i]]]$rd)
    }

    if(is.null(xmax)){
      xmax=rdo2@binParams$gainThresh*2 ##max(bins,na.rm=TRUE) ##bins[round(length(bins)*0.999)]
    }
    len=length(bins)
    e=rdo@entrypoints
    breaks=seq(0,xmax,windSize)
    oDisp=rdo@params$overDispersion
    
    med=rdo@binParams$med
    
    p1 = rpois.od(len*rdo@binParams$hapPerc,med/2,oDisp)
    p2=rpois.od(len*rdo@binParams$dipPerc,med,oDisp)
    p3 = rpois.od(len*rdo@binParams$tripPerc,(med*(3/2)),oDisp)
    
    p=append(append(p1,p2),p3)
    
    bins = bins[which(bins<xmax & bins>0)]
    p = p[which(p<xmax & p>0)]    

    plotOverlappingHist(a=bins,b=p,breaks=breaks)
    abline(v=rdo@binParams$gainThresh,col="red")
    abline(v=rdo@binParams$lossThresh,col="red")
  }

  dev.off()
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

writeSegs <- function(segs,rdo){
  write.table(segs, file=paste(rdo@params$outputDirectory,"/segs.dat",sep=""), sep="\t", quote=F, row.names=F, col.names=F)
}

writeAlts <- function(segs,rdo){
  write.table(getAlts(segs,rdo), file=paste(rdo@params$outputDirectory,"/alts.dat",sep=""), sep="\t", quote=F, row.names=F, col.names=F)
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



sort.data.frame <- function(x, by){
    # Author: Kevin Wright
    # with some ideas from Andy Liaw
    # http://tolstoy.newcastle.edu.au/R/help/04/07/1076.html
 
    # x: A data.frame
    # by: A one-sided formula using + for ascending and - for descending
    #     Sorting is left to right in the formula
  
    # Useage is:
    # library(nlme);
    # data(Oats)
    # sort(Oats, by= ~nitro-Variety)
 
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




calculateMedianFromDbsnpSites <- function(rdo, snpBinSize, peakWiggle=3, sites=NULL, plot=FALSE){
  ## read in dbsnp het sites
  if(is.null(sites)){
    if(verbose){
      print("reading in dbSNP het sites")
    }
    sites = read.table(rdo@params$dbSnpVaf)
  }

  sites = sites[which(((sites$V5 + sites$V6) <= 100) &
    ((sites$V5 + sites$V6) >= 20) &
    (sites$V7 > 15)),c(1,2,7)]

  homsites = sites[which(sites$V7 > 85),]
  hetsites = sites[which(sites$V7 < 85),]  


  ##free up some memory
  sites = NULL
  gc()

  
  if(verbose){
    print("finding clean windows of CN 2x, 3x, 4x")
  }
  
  
  doCalc <- function(rdo,snpBinSize,peakWiggle, chr, hetsites, homsites, plot){
    chrLen = getChrLength(chr,rdo@entrypoints)
    ## create an IRange for bins
    bins <- IRanges(start = (0:(ceiling(chrLen/snpBinSize)-1)*snpBinSize)+1, end = (1:ceiling(chrLen/snpBinSize))*snpBinSize)
    end(bins[length(bins)]) <- chrLen
        
    ##create IRange for dbsnp vafs
    homIsites = IRanges(start=homsites$V2, end=homsites$V2)
    hetIsites = IRanges(start=hetsites$V2, end=hetsites$V2)

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
      pdf(file=paste("vafplots/vafs.",chr,".pdf",sep=""))
    }


    for(i in 1:length(bins)){

      #if we have sites to look at
      if(length(hetsites[which(hetBinnedSites[,1]==i),]$V7) > 0){

        ##exclude sites that look 1x (low het to homo ratio)
        hetcount = length(hetsites[which(hetBinnedSites[,1]==i),1])
        homcount = length(homsites[which(homBinnedSites[,1]==i),1])

        if(hetcount/homcount < 0.75){          
          cn = NULL
        } else {
          den=(density(hetsites[which(hetBinnedSites[,1]==i),]$V7, bw=3))
          peaks=den$x[findPeaks(den$y)]
          peakHeight=den$y[findPeaks(den$y)]

          #filter out low, probably false peaks
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
      
    if(plot){
      dev.off()
    }
    
    if(verbose){      
      cat(paste("chr",chr," VAF anaylsis: \n  ",length(which(ploidy == 2))," 2x sites\n  ",length(which(ploidy == 3))," 3x sites\n  ",length(which(ploidy == 4))," 4x sites\n  ",length(bins)-length(ploidy)," sites ambiguous\n",sep=""))
    }
    
    ##now grab all the reads and match them up, adjust readcounts by ploidy, put the readcounts into a big list
    validWinds = IRanges(start=sts, end=sps)
    ##get all the positions and readcounts for this chr
    
    df = makedf(rdo@chrs,rdo@params)
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
    #print("test1")
    print(paste("mean:   ",mean(aReadDepths, na.rm=TRUE)))
    print(paste("median: ",median(aReadDepths, na.rm=TRUE)))
    return(aReadDepths)
  }
    

  
  ## find locations of 1MB windows that are
  ## unambiguously 2x, 3x, 4x

  ##this uses a ton of memory, so we're just going to use one core for now
  options(cores = 1);
  mcoptions <- list(preschedule = FALSE)
  adjReadDepths = foreach(chr=rdo@entrypoints$chr, .combine="append",.options.multicore=mcoptions) %dopar% {
    doCalc(rdo,snpBinSize,peakWiggle, chr, hetsites[which(hetsites$V1==chr),], homsites[which(homsites$V1==chr),], plot=plot)
  }
  ##restore cores for future use
  options(cores = rdo@params$maxCores)
  
  #return the average of the adj readcounts
  #print(adjReadDepths)
  #print("test2")
  #print(adjReadDepths)
  print(mean(adjReadDepths, na.rm=TRUE))
  pdf(file="vafplots/means.pdf")
  hist(adjReadDepths,breaks=100,col="darkgreen")
  abline(v=mean(adjReadDepths,na.rm=T),col="red")
  abline(v=median(adjReadDepths,na.rm=T),col="blue")
  dev.off()

  #return(adjReadDepths)
  return(mean(adjReadDepths, na.rm=TRUE));
}
