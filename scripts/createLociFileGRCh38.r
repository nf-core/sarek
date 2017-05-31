# cleaning the workspace
rm(list = ls(all = TRUE))

#Set the path to output folder to something relevent here:
setwd("/path/to/locifiles")

source("https://bioconductor.org/biocLite.R")
biocLite("biomaRt")
library(biomaRt)

variation = useEnsembl(biomart="ENSEMBL_MART_SNP", dataset="hsapiens_snp")
#listFIlters(variation)
ensembl_snps_GRCh38 <- getBM(attributes=c('chr_name','chrom_start','chrom_end','minor_allele_freq', 'consequence_type_tv'),
                      filters = "minor_allele_freq_second", 
                      values = "0.3", 
                      mart = variation)
#consequences=unique(ensembl_snps_GRCh38[,5])
loci_wgs_GRCh38=ensembl_snps_GRCh38[ensembl_snps_GRCh38$chrom_start == ensembl_snps_GRCh38$chrom_end,1:4]
write.table(loci_wgs_GRCh38, file="ensembl_GRCh38_maf0.3.loci", sep="\t", quote=F, row.names=F, col.names=F)
