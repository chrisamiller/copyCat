The copyCat package for R can detect somatic copy number aberrations by measuring the depth of coverage obtained by massively parallel sequencing of the genome. It achiev0;95;ces higher accuracy than many other packages, and runs faster by utilizing multi-core architectures to parallelize the processing of these large data sets.

copyCat takes in paired samples (tumor and normal) and can utilize mutation frequency information from samtools to help correct for purity and ploidy. This package also includes a method for effectively increasing the resolution obtained from low-coverage experiments by utilizing breakpoint information from paired end sequencing to do positional refinement.  It's primary input comes from running bam-window (https://github.com/genome-vendor/bam-window) on the tumor and normal bam files. 

#Installation

    #install devtools if you don't have it already
    install.packages("devtools")
    library(devtools)
    install_github("chrisamiller/copycat")

#Usage
    library(copyCat)
    #The most convenient way to run copyCat is through the functions in meta.R. 
    #For a paired tumor/normal sample, this looks something like this:
    runPairedSampleAnalysis(annotationDirectory="~/annotations/copyCat/hg19/",
                        outputDirectory="ccout",
                        normal="/path/to/normal_window_file
                        tumor="/path/to/tumor_window_file
                        inputType="bins",
                        maxCores=2,
                        binSize=0, #infer automatically from bam-window output
                        perLibrary=1, #correct each library independently
                        perReadLength=1, #correct each read-length independently
                        verbose=TRUE,
                        minWidth=3, #minimum number of consecutive winds need to call CN
                        minMapability=0.6, #a good default
                        dumpBins=TRUE,
                        doGcCorrection=TRUE,
                        samtoolsFileFormat="unknown", #will infer automatically - mpileup 10col or VCF
                        purity=1,
                        normalSamtoolsFile="normal_mpileup",
                        tumorSamtoolsFile="tumor_mpileup")  #uses the VAFs of mpileup SNPs to infer copy-neutral regions


#Annotations
CopyCat requires mapability and gc-content information that is dependent on the read-lengths of your data. (It accepts +/- 10bp as reasonable approximations) Annotation files that cover common read lengths on build37 are hosted here:
- [Base directory with genome info and anno for read length 100bp](https://dl.dropboxusercontent.com/u/21436449/copycat.anno.hg19.tar.gz)
- [read length 53](https://dl.dropboxusercontent.com/u/21436449/readlength.53.tar.gz)
- [read length 66](https://dl.dropboxusercontent.com/u/21436449/readlength.66.tar.gz)
- [read length 75](https://dl.dropboxusercontent.com/u/21436449/readlength.75.tar.gz)
- [read length 85](https://dl.dropboxusercontent.com/u/21436449/readlength.85.tar.gz)

#Notes
- The copyCat package is loosely based on [readDepth](https://code.google.com/p/readdepth/), a tool by the same author. 
- It does support single-sample CN calling using the "runSingleSampleAnalysis" function.
- It is not specific to the human genome. To create your own annotation files, use the above as a template, and fill in your own annoations for mapability (using self-aligments with your aligner of choice) and GC-content (for reads starting in each 100bp window).
- a window size of 10k is a reasonable default that balances specificity and sensitivity
