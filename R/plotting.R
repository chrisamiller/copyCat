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
##---------------------------------------------------------
## Internal function to wrap plotting
##
doSegPlot <- function(rdo,segs,chr){
  st = 1
  sp = rdo@entrypoints[which(rdo@entrypoints$chr == chr),]$length
  ymax=3 #multiplied by the median to get height
  med = rdo@binParams$med
  
  winds = rdo@chrs[[chr]]$rd
  binSize = rdo@binParams$binSize  
  pos = seq(binSize/2,(((length(winds)-1)*binSize)-binSize/2),binSize)
  pos = append(pos,sp)
  
  par(mar=c(5, 4, 4, 4) + 0.1)
  plot(pos,winds,ylab="number of reads", xlab="position (bp)", ylim=c(0,(median(winds,na.rm=T)*ymax)), pch=18, col=rgb(0,0,0,0.5), main=paste("Chr",chr,sep=""))

  abline(h=med,col="blue")
  abline(h=(rdo@binParams$gainThresh/2)*med,col=rgb(0,0.5,0,0.5))
  abline(h=(rdo@binParams$lossThresh/2)*med,col=rgb(0,0.5,0,0.5))

  asegs = segs[which(segs$chrom == chr),]
  for(i in 1:length(asegs[,1])){
    segments(asegs[i,2],(asegs[i,5]*(med/2)),asegs[i,3],(asegs[i,5]*(med/2)), col="red",lwd=3)
  }
  
  par(new=T)
  plot(-10000,-10000, ylim=c(0,(median(winds,na.rm=T)*(ymax*2))/med), xlim=c(1,sp), axes=F,xlab="", ylab="")
  axis(4, ylim=c(0,(median(winds,na.rm=T)*(ymax*2))/med), col="red",col.axis="red")
  abline(h=seq(0,100,1),col="grey50",lty=3)
  abline(v=seq(1,sp,2000000),col="grey50",lty=2)
  mtext("Copy Number",side=4,col="red",line=2.5)  

}

##----------------------------------------------
## plot the segments for a given chromosome
##
plotSegs <- function(rdo,segs,chr){
  pdf(paste(rdo@params$outputDirectory,"/plots/points.",chr,".pdf",sep=""),width=12,height=4)
  st = 1
  sp = rdo@entrypoints[which(rdo@entrypoints$chr == chr),]$length
  ymax=3 #multiplied by the median to get height
  med = rdo@binParams$med
  
  winds = rdo@chrs[[chr]]$rd
  binSize = rdo@binParams$binSize  
  pos = seq(binSize/2,(((length(winds)-1)*binSize)-binSize/2),binSize)
  pos = append(pos,sp)
  
  par(mar=c(5, 4, 4, 4) + 0.1)
  plot(pos,winds,ylab="number of reads", xlab="position (bp)", ylim=c(0,(median(winds,na.rm=T)*ymax)), pch=18, col=rgb(0,0,0,0.5), main=paste("Chr",chr,sep=""))

  abline(h=med,col="blue")
  abline(h=(rdo@binParams$gainThresh/2)*med,col=rgb(0,0.5,0,0.5))
  abline(h=(rdo@binParams$lossThresh/2)*med,col=rgb(0,0.5,0,0.5))

  asegs = segs[which(segs$chrom == chr),]
  for(i in 1:length(asegs[,1])){
    segments(asegs[i,2],(asegs[i,5]*(med/2)),asegs[i,3],(asegs[i,5]*(med/2)), col="red",lwd=3)
  }
  
  par(new=T)
  plot(-10000,-10000, ylim=c(0,(median(winds,na.rm=T)*(ymax*2))/med), xlim=c(1,sp), axes=F,xlab="", ylab="")
  axis(4, ylim=c(0,(median(winds,na.rm=T)*(ymax*2))/med), col="red",col.axis="red")
  abline(h=seq(0,100,1),col="grey50",lty=3)
  abline(v=seq(1,sp,2000000),col="grey50",lty=2)
  mtext("Copy Number",side=4,col="red",line=2.5)
  dev.off()
}


##----------------------------------------------
## plot the segments for a given chromosome
##
plotSegsPairedLog <- function(rdo.ref, rdo.test, segs,chr){
  pdf(paste(rdo.ref@params$outputDirectory,"/plots/points.",chr,".paired.log.pdf",sep=""),width=12,height=4)
  st = 1
  sp = rdo.ref@entrypoints[which(rdo.ref@entrypoints$chr == chr),]$length
  ymax=

  df = makeDfLogPaired(rdo.ref,rdo.test)
  df = df[df$chr==chr,]
  df$score = (2^df$score)*2
  par(mar=c(5, 4, 4, 4) + 0.1)
  plot(df$pos,df$score,ylab="log2 ratio", xlab="position (bp)", ylim=c(0,ymax), pch=18, col=rgb(0,0,0,0.5), main=paste("Chr",chr,sep=""),plot.first=abline(h=seq(0,100,1),col="grey50",lty=3))

  abline(h=2,col="blue")
  abline(h=rdo.ref@binParams$gainThresh,col="green")
  abline(h=rdo.ref@binParams$lossThresh,col="green")

  asegs = segs[which(segs$chrom == chr),]
  for(i in 1:length(asegs[,1])){
    segments(asegs[i,2], asegs[i,5], asegs[i,3], asegs[i,5], col="red", lwd=3)
  }
  
  abline(v=seq(1,sp,2000000),col="grey50",lty=2)
  
  dev.off()
}

##----------------------------------------------
## plot the segments for a given chromosome
##
plotSegsPaired <- function(rdo.ref, rdo.test, segs,chr){
  pdf(paste(rdo.ref@params$outputDirectory,"/plots/points.",chr,".paired.pdf",sep=""),width=12,height=8)
  
  st = 1
  sp = rdo.ref@entrypoints[which(rdo.ref@entrypoints$chr == chr),]$length
  ymax=3 #multiplier to median

  df = makeDfLogPaired(rdo.ref,rdo.test)
  df = df[df$chr==chr,]
  df$score = (2^df$score)*2
  par(mar=c(5, 4, 4, 4) + 0.1)
  plot(df$pos,df$score,ylab="log2 ratio", xlab="position (bp)", ylim=c(0,ymax), pch=18, col=rgb(0,0,0,0.5), main=paste("Chr",chr,sep=""),plot.first=abline(h=seq(0,100,1),col="grey50",lty=3))

  abline(h=2,col="blue")
  abline(h=rdo.ref@binParams$gainThresh,col="green")
  abline(h=rdo.ref@binParams$lossThresh,col="green")

  asegs = segs[which(segs$chrom == chr),]
  for(i in 1:length(asegs[,1])){
    segments(asegs[i,2], asegs[i,5], asegs[i,3], asegs[i,5], col="red", lwd=3)
  }
  
  abline(v=seq(1,sp,2000000),col="grey50",lty=2)
  
  dev.off()
}
