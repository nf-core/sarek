#!/bin/bash
set -xeuo pipefail

function nf_test() {
  nextflow run "$@" -profile travis
}

nf_test buildReferences.nf --download

# Clean up docker images
docker rmi -f maxulysse/igvtools:1.1
nf_test . --test --step preprocessing --tools MultiQC
# Clean up docker images
docker rmi -f maxulysse/fastqc:1.1 maxulysse/mapreads:1.1 maxulysse/picard:1.1
nf_test . --step realign
nf_test . --step realign --tools HaplotypeCaller -resume
nf_test . --step recalibrate
nf_test . --step skipPreprocessing --tools FreeBayes,HaplotypeCaller,MultiQC,MuTect1,MuTect2,Strelka
# Clean up docker images
docker rmi -f maxulysse/concatvcf:1.1 maxulysse/freebayes:1.1 maxulysse/gatk:1.1 maxulysse/mutect1:1.1 maxulysse/samtools:1.1 maxulysse/strelka:1.1
nf_test . --step skipPreprocessing --tools MuTect2,snpEff -resume
nf_test . --step annotate --tools MultiQC,snpEff --annotateTools MuTect2
nf_test . --step annotate --tools snpEff --annotateVCF VariantCalling/MuTect2/mutect2_9876T_vs_1234N.vcf.gz,VariantCalling/MuTect2/mutect2_9877R_vs_1234N.vcf.gz
