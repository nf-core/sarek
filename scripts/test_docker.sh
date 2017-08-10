#!/bin/bash
set -xeuo pipefail

function nf_test() {
  echo "$(tput setaf 1)nextflow run $@ -profile travis -resume --verbose$(tput sgr0)"
  nextflow run "$@" -profile travis -resume --verbose
}

nf_test buildReferences.nf --download

# Clean up docker images
docker rmi -f maxulysse/igvtools:1.1
nf_test . --test --step preprocessing
# Clean up docker images
docker rmi -f maxulysse/fastqc:1.1 maxulysse/mapreads:1.1 maxulysse/picard:1.1
nf_test . --step realign --noReports
nf_test . --step realign --tools HaplotypeCaller
nf_test . --step realign --tools HaplotypeCaller --noReports --noGVCF
nf_test . --step recalibrate --noReports
nf_test . --step recalibrate --tools FreeBayes,HaplotypeCaller,MuTect1,MuTect2,Strelka
# Test whether restarting from an already recalibrated BAM works
nf_test . --step skipPreprocessing --tools Strelka --noReports
# Clean up docker images
docker rmi -f maxulysse/concatvcf:1.1 maxulysse/freebayes:1.1 maxulysse/gatk:1.1 maxulysse/mutect1:1.1 maxulysse/samtools:1.1 maxulysse/strelka:1.1
nf_test . --step skipPreprocessing --tools MuTect2,snpEff,VEP --noReports
nf_test . --step annotate --tools snpEff,VEP --annotateTools MuTect2
nf_test . --step annotate --tools snpEff,VEP --annotateVCF VariantCalling/MuTect2/mutect2_9876T_vs_1234N.vcf.gz,VariantCalling/MuTect2/mutect2_9877R_vs_1234N.vcf.gz --noReports
