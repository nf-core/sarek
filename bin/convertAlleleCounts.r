#!/bin/env Rscript
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
################################################################################

##First read in the arguments listed at the command line
args = commandArgs(trailingOnly=TRUE)

## args is now a list of character vectors
## First check to see if arguments are passed.
if(length(args)<5){
    stop("No input files supplied\n\nUsage:\nRscript convertAlleleCounts.r tumorid tumorac normalid normalac gender\nWhere:\ntumorid - id of tumor sample\ntumorac - output from AlleleCount for the tumor\nnormalid - id of normal sample\nnormalac - output from AlleleCount for the normal\ngender - XX or XY\n\n")
} else{
    tumorid = args[1]
    tumorac = args[2]
    normalid = args[3]
    normalac = args[4]
    gender = args[5]
}

tumorcounts = read.table(tumorac, header=F, sep="\t")
normalcounts = read.table(normalac, header=F, sep="\t")

SNPpos = matrix(nrow = dim(normalcounts)[1],ncol = 2)

rownames(SNPpos) = paste("snp",1:dim(SNPpos)[1],sep="")

#Change rownames to "chr_pos" instead, such as 1_44552
#This does not exactly work:
#rownames(SNPpos) = apply(cbind(tumorcounts[,1], tumorcounts[,2]), 1, paste, collapse="_")
#This is for compatibility with gc correction file

colnames(SNPpos) = c("Chr","Position")
SNPpos[,1] = as.vector(normalcounts[,1])
SNPpos[,2] = normalcounts[,2]

#Caclulate BAF
Tumor_BAF = matrix(nrow = dim(normalcounts)[1],ncol = 1)
rownames(Tumor_BAF) = rownames(SNPpos)
colnames(Tumor_BAF) = c(tumorid)
acgt = tumorcounts[,c(3:6)]
acgts = t(apply(acgt,1,sort))
Tumor_BAF[,1] = acgts[,4]/(acgts[,3]+acgts[,4])
Tumor_BAF[,1] = ifelse(runif(length(Tumor_BAF[,1]))<0.5,Tumor_BAF[,1],1-Tumor_BAF[,1])
Tumor_BAF[is.nan(Tumor_BAF)]=NA

Germline_BAF = matrix(nrow = dim(normalcounts)[1],ncol = 1)
rownames(Germline_BAF) = rownames(SNPpos)
colnames(Germline_BAF) = c(normalid)
acgt = normalcounts[,c(3:6)]
acgts = t(apply(acgt,1,sort))
Germline_BAF[,1] = acgts[,4]/(acgts[,3]+acgts[,4])
Germline_BAF[,1] = ifelse(runif(length(Germline_BAF[,1]))<0.5,Germline_BAF[,1],1-Germline_BAF[,1])
Germline_BAF[is.nan(Germline_BAF)]=NA


Tumor_LogR = matrix(nrow = dim(normalcounts)[1],ncol = 1)
Germline_LogR = matrix(nrow = dim(normalcounts)[1],ncol = 1)
rownames(Tumor_LogR) = rownames(SNPpos)
colnames(Tumor_LogR) = c(tumorid)
rownames(Germline_LogR) = rownames(SNPpos)
colnames(Germline_LogR) = c(normalid)
Tumor_LogR[,1] = log(tumorcounts[,7]/normalcounts[,7],2)
Germline_LogR[,1] = 0
Tumor_LogR[is.infinite(Tumor_LogR)]=NA
if(gender=="XY") {
    Tumor_LogR[SNPpos[,1]=="X",1] = Tumor_LogR[SNPpos[,1]=="X",1]-1
    Germline_LogR[SNPpos[,1]=="X",1] = Germline_LogR[SNPpos[,1]=="X",1]-1
}
Tumor_LogR[,1] = Tumor_LogR[,1] - median(Tumor_LogR[,1],na.rm=T)
# set regions with 0 reads in tumor and normal to a LogR of 0.
Tumor_LogR[is.na(Tumor_LogR[,1]),1] = 0

# limit the number of digits:
Tumor_LogR = round(Tumor_LogR,4)
Tumor_BAF = round(Tumor_BAF,4)
Germline_LogR = round(Germline_LogR,4)
Germline_BAF = round(Germline_BAF,4)

# write output to files
write.table(cbind(SNPpos,Tumor_LogR),paste(tumorid,".LogR",sep=""),sep="\t",row.names=T,col.names=NA,quote=F)
write.table(cbind(SNPpos,Tumor_BAF),paste(tumorid,".BAF",sep=""),sep="\t",row.names=T,col.names=NA,quote=F)
write.table(cbind(SNPpos,Germline_LogR),paste(normalid,".LogR",sep=""),sep="\t",row.names=T,col.names=NA,quote=F)
write.table(cbind(SNPpos,Germline_BAF),paste(normalid,".BAF",sep=""),sep="\t",row.names=T,col.names=NA,quote=F)
