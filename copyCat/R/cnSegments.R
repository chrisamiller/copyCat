##----------------------------------------------------
## Use circular binary segmentation to merge adjacent
## windows and identify breakpoints
##
cnSegments <- function(rdo,onlyAlts=FALSE,minWidth=3,alpha=0.01,undoSD=2,rmGaps=TRUE){
  library('DNAcopy')
  df = makeDfLog(rdo@chrs,rdo@binParams)
  return(getSegs(rdo, df,rdo@binParams,rdo@entrypoints,onlyAlts,minWidth,alpha,rmGaps,undoSD))
}

##----------------------------------------------------
## Use circular binary segmentation to merge adjacent
## windows and identify breakpoints from tumor/normal samples
##
cnSegments.paired <- function(rdo.ref,rdo.test,onlyAlts=FALSE,minWidth=3,alpha=0.01,undoSD=2,rmGaps=TRUE){
  library('DNAcopy')

  if(verbose){
    ## print("downsampling genome with greater coverage, combining data to get log2 ratios")
    print("combining data to get log2 ratios")
  }

  df = makeDfLogPaired(rdo.ref,rdo.test)

  #do the segmentation
  return(getSegs(rdo.ref,df,params=rdo.ref@binParams, entrypoints=rdo.ref@entrypoints, onlyAlts=onlyAlts, minWidth=minWidth,
                 alpha=alpha, rmGaps=rmGaps, undoSD=undoSD))
}


##----------------------------------------------
## run circular binary segmentation to identify discrete
## segments of gain and loss
##
getSegs <- function(rdo, gd2, params, entrypoints, onlyAlts, minWidth=3, alpha=0.01, rmGaps=TRUE, undoSD=2){
  CNA.object <-CNA( genomdat = gd2$score, chrom = gd2$chr, maploc =gd2$pos, data.type = 'logratio')

  ## we could enable smoothing, but it seems to work better
  ## without this step.  Our data is less noisy than microarrays
  ## smoothed.CNA.object <-smooth.CNA(CNA.object)

  if(verbose){
    print("Segmenting Data");
  }
  #run cbs
  segs <- segment(CNA.object, verbose=0, alpha=alpha, min.width=minWidth, undo.splits="sdundo",undo.SD=undoSD)
  segs = segs$output
  if(rmGaps==TRUE){
    ##adjust first and last segs to meet chromosome ends,
    ##the remove gaps between segments by setting the breakpoint
    ##halfway between the calls.

    ## do not remove gaps if they are larger that this value
    ## (prevents creation of huge segments across centromeres)
    maxDist=5000000

    ##function to do the removal
    removeGaps <- function(chr,asegs,entrypoints){
      asegs = asegs[which(asegs$chrom == chr),]
      if(length(asegs[,1])==0){
        return(NULL)
      }
      entrypoint = entrypoints[which(entrypoints$chr==chr),]
      asegs[1,]$loc.start = 1
      asegs[length(asegs[,1]),]$loc.end = entrypoint$length
      if(length(asegs[,1])==1){
        return(asegs)
      }

      for(i in 1:(length(asegs[,1])-1)){
        ##don't merge if they're really far apart
        if((asegs[i+1,]$loc.start-asegs[i,]$loc.end) < maxDist){
          half = floor((asegs[i,]$loc.end + asegs[i+1,]$loc.start)/2)
          asegs[i,]$loc.end = half
          asegs[i+1,]$loc.start = half + 1
        } else {
          if(verbose){
            cat("not removing",chr,asegs[i,]$loc.end,asegs[i+1,]$loc.start," - gap too large\n")
          }
        }
      }
      return(asegs)
    }
    ##do the gap removal
    segs2 = NULL
    for(i in 1:length(entrypoints$chr)){
      chr = entrypoints$chr[i]
      chrSegs = removeGaps(chr,segs,entrypoints)
      if(!is.null(chrSegs)){
        if(is.null(segs2)){
          segs2 = chrSegs
        }else{
          segs2 = rbind(segs2,chrSegs)
        }
      }
    }
    segs = segs2
  }


  ##convert means back out of log space into absolute copy number
  segs$seg.mean=(2^segs$seg.mean)*2

##   ## alternate plots from the DNAcopy package
##   ## disabled for now
##   doPlots <- function(){
##     pdf("output/CNplots.pdf")
##     plot(segs,plot.type="w")
##     plot(segs,plot.type="s")
##     plot(segs,plot.type="p")

##     for(i in unique(segs$output$chrom)){
##       a = subset(segs,chromlist=c(i))
##       plot(a)
##       mtext(paste("Chromosome",i))
##     }
##     dev.off()
##   }
##   doPlots()

  if(onlyAlts){
    return(getAlts(segs[,2:6],rdo))
  }
  ##remove the sample name column before returning
  return(segs[,2:6])
}


##----------------------------------------------
## merge segments if they are adjacent and both gains or losses
## useful for some types of plotting
mergeSegs <- function(segs,rdo){
  df = NULL
  params=rdo@binParams
  gainCN = params$gainThresh/params$med
  lossCN = params$lossThresh/params$med

  buffer = segs[1,]
  for(i in 2:length(segs[,1])){
    ##if they are adjacent
    if(buffer[3]+1 == segs[i,2]){
      ## if they're both gains or losses
      if(((buffer[5] > gainCN) & (segs[i,5] > gainCN)) |
         ((buffer[5] > lossCN) & (segs[i,5] > lossCN))){

        #merge them
        cat("merge\n")
        buffer[3] = segs[i,3]
        buffer[5] = ((buffer[4]*buffer[5])+(segs[i,4]*segs[i,5]))/(buffer[4]+segs[i,4])
        buffer[4] = buffer[4] + segs[i,4]

      }else{
        df = rbind(df,buffer)
        buffer = segs[i,]
      }
    }else{
      df = rbind(df,buffer)
      buffer = segs[i,]
    }
  }
  #clear buffer at end
  df = rbind(df,buffer)
  return(df)
}



##-----------------------------------------------
## subset out just the alterations from the copy number segments
##
getAlts <- function(segs,rdo){
  return(segs[which(segs$seg.mean > rdo@binParams$gainThresh | segs$seg.mean < rdo@binParams$lossThresh),])
}


##--------------------------------------------------
## trim the ends of segments such that they're called
## normal (2.0) until the first bin that isn't NA
##
## useful for chromosomes like 22, where there are tens of
## megabases that are unmappable at beginning. Without this,
## those regions would be called with the same CN as the first
## few mapable bins
##
trimSegmentEnds <- function(segs,rdo){

  doTrimming <- function(chr, asegs, bins){
    ##find first pos that doesn't have an NA
    st=1
    while(is.na(bins[st,"rd"]) & (st < length(bins[,1]))){
      st = st + 1
    }

    ##find last pos that doesn't have an NA
    sp=length(bins[,1])
    while(is.na(bins[sp,"rd"]) & (sp > 0)){
      sp = sp - 1
    }

    ##edge case 1 - everything is NA, do nothing
    if(sp == 0){
      return(asegs)
    }

    ##trim start (unless there are no NAs at beginning)
    if(!(st == 1)){
      asegs[1,]$loc.start=(st-1)*rdo@binParams$binSize+1
      a=data.frame(chrom=chr,
        loc.start=1,
        loc.end=(st-1)*rdo@binParams$binSize,
        num.mark=0,
        seg.mean=2.0)

      asegs=rbind(a,asegs)
    }

    ##trim end (unless there are no NAs at end)
    if(!(sp == length(bins[,1]))){
      asegs[length(asegs[,1]),]$loc.end=(sp-1)*rdo@binParams$binSize-1
      a=data.frame(chrom=chr,
        loc.start=sp*rdo@binParams$binSize,
        loc.end=rdo@entrypoints[which(rdo@entrypoints$chr==chr),]$length,
        num.mark=0,
        seg.mean=2.0)
        asegs=rbind(asegs,a)
    }
    return(asegs)
  }

  chr = NULL
  foreach(chr=rdo@entrypoints$chr, .combine="rbind") %dopar% {
    doTrimming(chr, segs[which(segs$chrom==chr),], rdo@chrs[[chr]])
  }
}
