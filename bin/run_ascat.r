#!//bin/env Rscript


#######################################################################################################
# Description:
# R-script for converting output from AlleleCount to BAF and LogR values.
#
# Input:
# AlleleCounter output file for tumor and normal samples
# The first line should contain a header describing the data
# The following columns and headers should be present:
# CHR    POS     Count_A Count_C Count_G Count_T Good_depth
#
# Output:
# BAF and LogR tables (tab delimited text files)
#######################################################################################################

source("./ascat.R")
#gcfile = "/proj/b2011185/nobackup/wabi/somatic_wgs_pipeline/1000G_maf0.3_gc_id.txt"

args = commandArgs(trailingOnly=TRUE)

##args is now a list of character vectors
## First check to see if arguments are passed.
if(length(args)==0){
    stop("No input files supplied\n\nUsage:\nRscript run_ascat.r tumor_baf tumor_logr normal_baf normal_logr\n\n")
} else{
    tumorbaf = args[1]
    tumorlogr = args[2]
    normalbaf = args[3]
    normallogr = args[4]
}

#Load the  data
ascat.bc <- ascat.loadData(Tumor_LogR_file=tumorlogr, Tumor_BAF_file=tumorbaf, Germline_LogR_file=normallogr, Germline_BAF_file=normalbaf)

#gc correction
#This doens not wurk unless .BAF and .LogR files have been filtered with script "filter_logr_baf.py" so that
#the row names are on format chr_pos instead of snp1, snp2 etc.
#ascat.bc = ascat.GCcorrect(ascat.bc, GCcontentfile=gcfile)

#Plot the raw data
ascat.plotRawData(ascat.bc)

#Segment the data
ascat.bc <- ascat.aspcf(ascat.bc)

#Plot the segmented data
ascat.plotSegmentedData(ascat.bc)

#Run ASCAT to fit every tumor to a model, inferring ploidy, normal cell contamination, and discrete copy numbers
ascat.output <- ascat.runAscat(ascat.bc)
str(ascat.output)
plot(sort(ascat.output$aberrantcellfraction))
plot(density(ascat.output$ploidy))


