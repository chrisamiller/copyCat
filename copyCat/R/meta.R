##-----------------------------------------------------------
## pull out the subset of reads with a specific read length
##
runSingleSampleAnalysis <- function(annotationDirectory, outputDirectory, inputFile,
                                    inputType, maxCores=0, #0 means use all available cores
                                    binSize=0,  #0 means let copyCat choose (or infer from bins file)
                                    gcWindowSize=100, fdr=0.01, perLibrary=TRUE,
                                    perReadLength=TRUE, readLength=0,
                                    samtoolsFile=NULL, peakwiggle=3, snpBinSize=1000000,
                                    verbose=TRUE){

  verbose <<- verbose
  ##create the object
  rdo = new("rdObject")
  ##set the parameters
  rdo = setParams(rdo, annotationDirectory=annotationDirectory,
    outputDirectory=outputDirectory, inputType=inputType, maxCores=maxCores,
    inputFile=inputFile, binSize=binSize, perReadLength=perReadLength,
    perLibrary=perLibrary, readLength=readLength, gcWindowSize=gcWindowSize,
    fdr=fdr, maxCores=maxCores)
  ##bin the reads
  rdo=getReadDepth(rdo)
  ##correct for mapability
  rdo=mapCorrect(rdo)
  ##correct for gc-content
  rdo=gcCorrect(rdo)
  ##merge the corrected counts into one column
  rdo=mergeLibraries(rdo)

  ##use the samtools file to get cn-neutral median read count
  if(!(is.null(samtoolsFile))){
    rdo@params$med = cnNeutralDepthFromHetSites(rdo,samtoolsFile,1000000,peakWiggle=3,plot=TRUE)
  } else {
    rdo@params@med = getMedianReadCount(rdo)
  }

  ##segment the data using CBS
  segs = cnSegments(rdo)
  ##set gain and loss thresholds  
  rdo@binParams$gainThresh = 2.5
  rdo@binParams$lossThresh = 1.5
  ##write some output
  writeSegs(segs,rdo)
  writeAlts(segs,rdo)
}



##-----------------------------------------------------------
## pull out the subset of reads with a specific read length
##
runPairedSampleAnalysis <- function(annotationDirectory, outputDirectory, normal, tumor,
                                    inputType, maxCores=0, #0 means use all available cores
                                    binSize=0,  #0 means let copyCat choose (or infer from bins file)
                                    gcWindowSize=100, fdr=0.01, perLibrary=TRUE,
                                    perReadLength=TRUE, readLength=0, verbose=TRUE,
                                    outputSingleSample=FALSE, tumorSamtoolsFile=NULL,
                                    dumpBins=FALSE, minWidth=3){

  verbose <<- verbose

  ##create the object
  rdo = new("rdObject")
  ##set the parameters
  rdo = setParams(rdo, annotationDirectory=annotationDirectory,
    outputDirectory=outputDirectory, inputType=inputType,
    inputFile=normal, binSize=binSize, perReadLength=perReadLength,
    perLibrary=perLibrary, readLength=readLength, maxCores=maxCores)
  rdo@params$prefix="normal";
  
  ##bin the reads
  rdo=getReadDepth(rdo)

  ##no correction for mapability, but still remove low-mapability regions
  ##from consideration 
  rdo=mapCorrect(rdo,skipCorrection=TRUE);  
  ##correct for gc-content
  rdo=gcCorrect(rdo)
  ##merge the corrected counts into one column
  rdo=mergeLibraries(rdo)

  ##create the object
  rdo2 = new("rdObject")
  ##set the parameters
  rdo2 = setParams(rdo2, annotationDirectory=annotationDirectory,
    outputDirectory=outputDirectory, inputType=inputType,
    inputFile=tumor, binSize=binSize, perReadLength=perReadLength,
    perLibrary=perLibrary, readLength=readLength, maxCores=maxCores)
  rdo2@params$prefix="tumor";

  ##bin the reads
  rdo2=getReadDepth(rdo2)

  ##no correction for mapability, but still remove low-mapability regions
  ##from consideration 
  rdo2=mapCorrect(rdo2,skipCorrection=TRUE);
  ##correct for gc-content
  rdo2=gcCorrect(rdo2)
  ##merge the corrected counts into one column
  rdo2=mergeLibraries(rdo2)

  if(dumpBins){
    writePairedBins(rdo,rdo2)
  }
  
  ##segment the paired data using CBS
  segs = cnSegments.paired(rdo, rdo2, minWidth=minWidth)
  ##set gain and loss thresholds
  rdo@binParams$gainThresh = 2.5
  rdo@binParams$lossThresh = 1.5
  ##write some output
  writeSegs(segs,rdo,"segs.paired.dat")
  writeAlts(segs,rdo,"alts.paired.dat")


  if(outputSingleSample){
    #tumor
    if(!(is.null(tumorSamtoolsFile))){
      rdo2@params$med = cnNeutralDepthFromHetSites(rdo2,tumorSamtoolsFile,1000000,peakWiggle=3,plot=TRUE)
    } else {
      rdo2@params@med = getMedianReadCount(rdo2)
    }
    ##segment the paired data using CBS
    segs = cnSegments(rdo2)
    ##set gain and loss thresholds
    rdo2@binParams$gainThresh = 2.5
    rdo2@binParams$lossThresh = 1.5
    ##write some output
    writeSegs(segs,rdo,"segs.singleSample.tumor.dat")
    writeAlts(segs,rdo,"alts.singleSample.tumor.dat")

    
    #normal
    if(!(is.null(tumorSamtoolsFile))){
      rdo@params$med = cnNeutralDepthFromHetSites(rdo,tumorSamtoolsFile,1000000,peakWiggle=3,plot=TRUE)
    } else {
      rdo@params@med = getMedianReadCount(rdo)
    }
    ##segment the paired data using CBS
    segs = cnSegments(rdo)
    ##set gain and loss thresholds
    rdo@binParams$gainThresh = 2.5
    rdo@binParams$lossThresh = 1.5
    ##write some output
    writeSegs(segs,rdo,"segs.singleSample.normal.dat")
    writeAlts(segs,rdo,"alts.singleSample.normal.dat")
    
  }
  
  dumpParams(rdo)
  dumpParams(rdo2)
    
}
