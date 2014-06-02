The copyCat package for R can detect somatic copy number aberrations by measuring the depth of coverage obtained by massively parallel sequencing of the genome. It achieves higher accuracy than many other packages, and runs faster by utilizing multi-core architectures to parallelize the processing of these large data sets.

copyCat currently requires paired samples (tumor and normal) and can utilize mutation frequency information from samtools to help correct for purity and ploidy. This package also includes a method for effectively increasing the resolution obtained from low-coverage experiments by utilizing breakpoint information from paired end sequencing to do positional refinement.  It's primary input comes from running bam-window (https://github.com/genome-vendor/bam-window) on the tumor and normal bam files. 

#Installation

    #install devtools if you don't have it already
    install.packages("devtools")
    library(devtools)
    install_github("chrisamiller/copycat")


copyCat is loosely based on readDepth, a previous package by the same author.

