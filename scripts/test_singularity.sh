#!/bin/bash
set -xeuo pipefail

function nf_test() {
  echo "$(tput setaf 1)nextflow run $@ -profile singularityTest$(tput sgr0)"
  nextflow run "$@" -profile singularityTest
}

nf_test buildReferences.nf --download

# Clean up docker images
nf_test . --test --step preprocessing --reports
# Clean up docker images
nf_test . --step realign
nf_test . --step realign --tools HaplotypeCaller --reports -resume
nf_test . --step recalibrate --tools FreeBayes,HaplotypeCaller,MuTect1,MuTect2,Strelka --reports
# Test whether restarting from an already recalibrated BAM works
nf_test . --step skipPreprocessing --tools Strelka
# Clean up docker images
nf_test . --step skipPreprocessing --tools MuTect2,snpEff -resume
nf_test . --step annotate --tools snpEff --annotateTools MuTect2 --reports
nf_test . --step annotate --tools snpEff --annotateVCF VariantCalling/MuTect2/mutect2_9876T_vs_1234N.vcf.gz,VariantCalling/MuTect2/mutect2_9877R_vs_1234N.vcf.gz
