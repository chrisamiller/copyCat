##-----------------------------------------------------------
## pull out the subset of reads with a specific read length
##
runSingleSampleAnalysis <- function(annotationDirectory, outputDirectory, inputFile,
                                    inputType, maxCores=0, #0 means use all available cores
                                    binSize=0,  #0 means let readDepth choose (or infer from bins file)
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
    fdr=fdr)
  ##bin the reads
  rdo=readDepth(rdo)
  ##correct for mapability
  rdo=rd.mapCorrect2(rdo)
  ##correct for gc-content
  rdo=rd.gcCorrect(rdo)
  ##merge the corrected counts into one column
  rdo=mergeLibraries(rdo)

  ##use the samtools file to get cn-neutral median read count
  if(!(is.null(samtoolsFile))){
    rdo@params$med = cnNeutralDepthFromHetSites(rdo,samtoolsFile,1000000,peakWiggle=3,plot=TRUE)
  } else {
    rdo@params@med = getMedianReadCount(rdo)
  }

  ##segment the data using CBS
  segs = rd.cnSegments(rdo)
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
                                    binSize=0,  #0 means let readDepth choose (or infer from bins file)
                                    gcWindowSize=100, fdr=0.01, perLibrary=TRUE,
                                    perReadLength=TRUE, readLength=0, verbose=TRUE,
                                    outputSingleSample=FALSE, nrmSamtoolsFile=NULL,
                                    tumorSamtoolsFile=NULL){

  verbose <<- verbose

  ##create the object
  rdo = new("rdObject")
  ##set the parameters
  rdo = setParams(rdo, annotationDirectory=annotationDirectory,
    outputDirectory=outputDirectory, inputType=inputType,
    inputFile=normal, binSize=binSize, perReadLength=perReadLength,
    perLibrary=perLibrary, readLength=readLength)
  rdo2@params$prefix="normal";
  
  ##bin the reads
  rdo=readDepth(rdo)
  ##correct for mapability
  rdo=rd.mapCorrect2(rdo)
  ##correct for gc-content
  rdo=rd.gcCorrect(rdo)
  ##merge the corrected counts into one column
  rdo=mergeLibraries(rdo)


  ##create the object
  rdo2 = new("rdObject")
  ##set the parameters
  rdo2 = setParams(rdo2, annotationDirectory=annotationDirectory,
    outputDirectory=outputDirectory, inputType=inputType,
    inputFile=tumor, binSize=binSize, perReadLength=perReadLength,
    perLibrary=perLibrary, readLength=readLength)
  rdo2@params$prefix="tumor";

  ##bin the reads
  rdo2=readDepth(rdo2)
  ##correct for mapability
  rdo2=rd.mapCorrect2(rdo2)
  ##correct for gc-content
  rdo2=rd.gcCorrect(rdo2)
  ##merge the corrected counts into one column
  rdo2=mergeLibraries(rdo2)


  ##segment the paired data using CBS
  segs = rd.paired.cnSegments(rdo,rdo2)
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
    segs = rd.cnSegments(rdo2)
    ##set gain and loss thresholds
    rdo2@binParams$gainThresh = 2.5
    rdo2@binParams$lossThresh = 1.5
    ##write some output
    writeSegs(segs,rdo,"segs.singleSample.tumor.dat")
    writeAlts(segs,rdo,"alts.singleSample.tumor.dat")

    
    #normal
    if(!(is.null(tumorSamtoolsFile))){
      rdo@params$med = cnNeutralDepthFromHetSites(rdo,normalSamtoolsFile,1000000,peakWiggle=3,plot=TRUE)
    } else {
      rdo@params@med = getMedianReadCount(rdo)
    }
    ##segment the paired data using CBS
    segs = rd.cnSegments(rdo)
    ##set gain and loss thresholds
    rdo@binParams$gainThresh = 2.5
    rdo@binParams$lossThresh = 1.5
    ##write some output
    writeSegs(segs,rdo,"segs.singleSample.normal.dat")
    writeAlts(segs,rdo,"alts.singleSample.normal.dat")
    
  }
  
  
}
