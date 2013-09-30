library(copyCat)
#do a paired sample analysis on the test data
runPairedSampleAnalysis(annotationDirectory="annotations/hg19.chr14only",
                        outputDirectory="results",
                        normal="data/normal.wind",
                        tumor="data/tumor.wind",
                        inputType="bins",
                        maxCores=4,
                        binSize=0,
                        perLibrary=1,
                        perReadLength=1,
                        verbose=TRUE,
                        dumpBins=TRUE,
                        doGcCorrection=TRUE,
                        normalSamtoolsFile="data/normal.samtools",
                        tumorSamtoolsFile="data/tumor.samtools",
                        purity=0.95)
print("test completed");
