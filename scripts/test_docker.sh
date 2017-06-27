#!/bin/bash
set -xeuo pipefail

function nf_test() {
  echo "$(tput setaf 1)nextflow run $@ -profile travis$(tput sgr0)"
  nextflow run "$@" -profile travis
}

nf_test buildReferences.nf --download

# Clean up docker images
docker rmi -f maxulysse/igvtools:1.1
nf_test . --test --step preprocessing --reports
# Clean up docker images
docker rmi -f maxulysse/fastqc:1.1 maxulysse/mapreads:1.1 maxulysse/picard:1.1
nf_test . --step realign
nf_test . --step realign --tools HaplotypeCaller --reports -resume
nf_test . --step recalibrate --tools FreeBayes,HaplotypeCaller,MuTect1,MuTect2,Strelka --reports
# Test whether restarting from an already recalibrated BAM works
nf_test . --step skipPreprocessing --tools Strelka
# Clean up docker images
docker rmi -f maxulysse/concatvcf:1.1 maxulysse/freebayes:1.1 maxulysse/gatk:1.1 maxulysse/mutect1:1.1 maxulysse/samtools:1.1 maxulysse/strelka:1.1
nf_test . --step skipPreprocessing --tools MuTect2,snpEff -resume
nf_test . --step annotate --tools snpEff --annotateTools MuTect2 --reports
nf_test . --step annotate --tools snpEff --annotateVCF VariantCalling/MuTect2/mutect2_9876T_vs_1234N.vcf.gz,VariantCalling/MuTect2/mutect2_9877R_vs_1234N.vcf.gz
