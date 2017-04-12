#!/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
if(length(args)<6){
    stop("No input files supplied\n\nUsage:\nRscript run_ascat.r tumor_baf tumor_logr normal_baf normal_logr tumor_sample_name baseDir\n\n")
} else{
    tumorbaf = args[1]
    tumorlogr = args[2]
    normalbaf = args[3]
    normallogr = args[4]
    tumorname = args[5]
    baseDir = args[6]

}

.libPaths((paste(baseDir,"/scripts", sep="") .libPaths()))
source(paste(baseDir,"/scripts/ascat.R", sep=""))

if(!require(RColorBrewer)){
    source("http://bioconductor.org/biocLite.R")
    biocLite("RColorBrewer", suppressUpdates=TRUE, lib="$baseDir/scripts")
    library(RColorBrewer)
}
options(bitmapType='cairo')


#Load the  data
ascat.bc <- ascat.loadData(Tumor_LogR_file=tumorlogr, Tumor_BAF_file=tumorbaf, Germline_LogR_file=normallogr, Germline_BAF_file=normalbaf)

#Plot the raw data
ascat.plotRawData(ascat.bc)

#Segment the data
ascat.bc <- ascat.aspcf(ascat.bc)

#Plot the segmented data
ascat.plotSegmentedData(ascat.bc)

#Run ASCAT to fit every tumor to a model, inferring ploidy, normal cell contamination, and discrete copy numbers
ascat.output <- ascat.runAscat(ascat.bc)

#Write out segmented regions
write.table(ascat.output$segments, file=paste(tumorname, "tumor.segments.txt", sep=""), sep="\t", quote=F, row.names=F)
