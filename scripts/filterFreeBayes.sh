#!/bin/bash
# Filtering somatic variations from FreeBayes VCF results
# $1 name of the VCF file (uncompressed)
# $2 ID of the normal sample
# $3 ID of the tumour sample
SAREKDIR=/home/szilva/dev/forkCAW

module load vcflib
NORMAL=$2
TUMOUR=$3
vcfsamplediff -s VT ${NORMAL} ${TUMOUR} $1 | awk '/^#/{print}!/^#/&&!/germline/'  > diffed.vcf
CHROM_LINE=(`grep -m 1 CHROM diffed.vcf`)
# calculate the 0-based indexes
T_INDEX=$((`echo ${CHROM_LINE[@]/${TUMOUR}//} | cut -d/ -f1 | wc -w | tr -d ' '`+1))
N_INDEX=$((`echo ${CHROM_LINE[@]/${NORMAL}//} | cut -d/ -f1 | wc -w | tr -d ' '`+1))
vcffilter -f "QUAL > 20" diffed.vcf| vcfflatten > flattened.vcf
awk -f ${SAREKDIR}/scripts/speedseq.filter.awk -v NORMAL_IDX=${N_INDEX} -v TUMOUR_IDX=${T_INDEX} flattened.vcf> ${1%.vcf}.filtered.vcf
# you can remove the flattened and the diffed files 
