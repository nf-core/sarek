#!/usr/bin/env Rscript
library("optparse")
option_list = list(
    make_option("--tumorbaf", type="character", default=NULL,
              help="tumor BAF file", metavar="character"),
              make_option("--tumorlogr", type="character", default=NULL,
              help="tumor LogR file", metavar="character"),
              make_option("--normalbaf", type="character", default=NULL,
              help="normal BAF file", metavar="character"),
              make_option("--normallogr", type="character", default=NULL,
              help="normal LogR file", metavar="character"),
              make_option("--tumorname", type="character", default=NULL,
              help="name of tumor sample file", metavar="character"),
              make_option("--basedir", type="character", default=NULL,
              help="main Sarek directory for sample", metavar="character"),
              make_option("--gcfile", type="character", default=NULL,
              help="GC correction file", metavar="character"),
              make_option("--gender", type="character", default=NULL,
              help="gender on format XX or XY", metavar="character"),
              make_option("--purity", type="double", default=NULL,
              help="override Ascat purity parameter (rho_manual) ", metavar="character"),
              make_option("--ploidy", type="double", default=NULL,
              help="override Ascat ploidy parameter (psi_manual)", metavar="character"))

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

if(is.null(opt$tumorbaf) || is.null(opt$tumorlogr) || is.null(opt$normalbaf) || is.null(opt$normallogr) || is.null(opt$tumorname) || is.null(opt$basedir) || is.null(opt$gcfile) || is.null(opt$gender))  {
    print_help(opt_parser)
    stop("At least one of the required arguments missing.", call.=FALSE)
    } 
    
library(ASCAT)

if(!require(RColorBrewer)){
    source("http://bioconductor.org/biocLite.R")
    biocLite("RColorBrewer", suppressUpdates=TRUE, lib="$opt$basedir/scripts")
    library(RColorBrewer)
}
options(bitmapType='cairo')


#Load the  data
#ascat.bc <- ascat.loadData(Tumor_LogR_file=opt$tumorlogr, Tumor_BAF_file=opt$tumorbaf, Germline_LogR_file=opt$normallogr, Germline_BAF_file=opt$normalbaf, chrs = c(1:22,"X","Y"), opt$gender = opt$gender, sexchromosomes = c("X", "Y"))

if(opt$gender=="XY"){
    ascat.bc = ascat.loadData(Tumor_LogR_file = opt$tumorlogr, Tumor_BAF_file = opt$tumorbaf, Germline_LogR_file = opt$normallogr, Germline_BAF_file = opt$normalbaf, gender = opt$gender, chrs = c(1:22,"X","Y"), sexchromosomes = c("X","Y"))
} else {
    ascat.bc = ascat.loadData(Tumor_LogR_file = opt$tumorlogr, Tumor_BAF_file = opt$tumorbaf, Germline_LogR_file = opt$normallogr, Germline_BAF_file = opt$normalbaf, gender = opt$gender, chrs = c(1:22,"X"), sexchromosomes = c("X"))

}


#GC wave correction
ascat.bc = ascat.GCcorrect(ascat.bc, opt$gcfile)

#Plot the raw data
ascat.plotRawData(ascat.bc)

#Segment the data
ascat.bc <- ascat.aspcf(ascat.bc)

#Plot the segmented data
ascat.plotSegmentedData(ascat.bc)

#Run ASCAT to fit every tumor to a model, inferring ploidy, normal cell contamination, and discrete copy numbers
#If psi and rho are manually set:
if (!is.null(opt$purity) && !is.null(opt$ploidy)){
  ascat.output <- ascat.runAscat(ascat.bc, gamma=1, rho_manual=opt$purity, psi_manual=opt$ploidy)
} else if(!is.null(opt$purity) && is.null(opt$ploidy)){
  ascat.output <- ascat.runAscat(ascat.bc, gamma=1, rho_manual=opt$purity)
} else if(!is.null(opt$ploidy) && is.null(opt$purity)){
  ascat.output <- ascat.runAscat(ascat.bc, gamma=1, psi_manual=opt$ploidy)
} else {
  ascat.output <- ascat.runAscat(ascat.bc, gamma=1)
}


#Write out segmented regions (including regions with one copy of each allele)
#write.table(ascat.output$segments, file=paste(opt$tumorname, ".segments.txt", sep=""), sep="\t", quote=F, row.names=F)

#Write out CNVs in bed format
cnvs=ascat.output$segments[2:6]
write.table(cnvs, file=paste(opt$tumorname,".cnvs.txt",sep=""), sep="\t", quote=F, row.names=F, col.names=T)

#Write out purity and ploidy info
summary <- tryCatch({
		matrix(c(ascat.output$aberrantcellfraction, ascat.output$ploidy), ncol=2, byrow=TRUE)}, error = function(err) {
			# error handler picks up where error was generated
			print(paste("Could not find optimal solution:  ",err))
			return(matrix(c(0,0),nrow=1,ncol=2,byrow = TRUE))
		}
)
colnames(summary) <- c("AberrantCellFraction","Ploidy")
write.table(summary, file=paste(opt$tumorname,".purityploidy.txt",sep=""), sep="\t", quote=F, row.names=F, col.names=T)
