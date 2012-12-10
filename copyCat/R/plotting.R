##--------------------------------------------------
## plot the loess correction for gc content
##
plotGcLoess <- function(med,gc,reads,fitted,fixed,fixedfit,rdo,libNum,readLength,libName){
  params=rdo@params
  if(is.null(params$prefix)){
    pdf(paste(params$outputDirectory,"/plots/gccontent.lib",libNum,".readLength",readLength,".pdf",sep=""),width=8,height=11)
  } else {
    pdf(paste(params$outputDirectory,"/plots/",params$prefix,".gccontent.lib",libNum,".readLength",readLength,".pdf",sep=""),width=8,height=11)
  }
  par(mfcol=c(2,1))

  ymax <- med*3
  ##plot original
  plot(gc,reads,ylim=c(0,ymax),ylab="mean # reads", main="GC content bias - raw Data", xlab="GC content %")
  mtext(paste("lib:",names(rdo@chrs[[1]])[libNum]," readLength:",readLength,sep=""))
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
plotMapLoess <- function(med,map,reads,fitted,fixed,fixedfit,rdo,libNum,readLength,libName){
  params=rdo@params
  if(is.null(params$prefix)){
    pdf(paste(params$outputDirectory,"/plots/mapability.lib",libNum,".readLength",readLength,".pdf",sep=""),width=8,height=11)
  } else {
    pdf(paste(params$outputDirectory,"/plots/",params$prefix,".mapability.lib",libNum,".readLength",readLength,".pdf",sep=""),width=8,height=11)
  }
  par(mfcol=c(2,1))

  ymax <- med*3
  ##plot original
  plot(map,reads,ylim=c(0,ymax),ylab="mean # reads", main="Mapability bias - raw Data", xlab="Mapability content %")
  mtext(paste("lib:",names(rdo@chrs[[1]])[libNum]," readLength:",readLength,sep=""))
  points(map,fitted,col="green",pch=20)
  abline(h=med,col="red")

  ## plot post-adjustment
  plot(map,fixed,ylim=c(0,ymax),ylab="mean # reads", main="Mapability bias - post LOESS correction", xlab="Mapability content %")
  points(map,fixedfit,col="green",pch=20)
  abline(h=med,col="red")
  
  dev.off()  
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
tripDist <- function(xmax=NULL,filename,windSize=20, rdo2, rdo3, rdo4){
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
