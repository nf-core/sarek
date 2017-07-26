#!/bin/bash
set -xeuo pipefail

function nf_test() {
  echo "$(tput setaf 1)nextflow run $@ -profile singularityTest$(tput sgr0)"
  nextflow run "$@" -profile singularityTest
}

nf_test buildReferences.nf --download

nf_test . --test
nf_test . --step realign --noReports
nf_test . --step realign --tools HaplotypeCaller -resume
nf_test . --step recalibrate --tools FreeBayes,HaplotypeCaller,MuTect2,Strelka
# Test whether restarting from an already recalibrated BAM works
nf_test . --step variantcalling --tools Strelka --noReports
nf_test . --step variantcalling --tools MuTect2,snpEff -resume  --noReports
nf_test . --step annotate --tools snpEff --annotateTools MuTect2
nf_test . --step annotate --tools snpEff --annotateVCF VariantCalling/MuTect2/mutect2_9876T_vs_1234N.vcf.gz,VariantCalling/MuTect2/mutect2_9877R_vs_1234N.vcf.gz  --noReports
