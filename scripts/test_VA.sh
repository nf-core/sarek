#!/bin/bash
set -xeuo pipefail
PROFILE=$1

function nf_test() {
  echo "$(tput setaf 1)nextflow run $@ -profile $PROFILE -resume --verbose$(tput sgr0)"
  nextflow run "$@" -profile $PROFILE -resume --verbose
}

nf_test buildReferences.nf --download

# Clean up images
if [ $PROFILE = travis ]
then
  docker rmi -f maxulysse/igvtools:1.1
else
  rm -rf work/singularity/igvtools-1.1.img
fi

nf_test . --step preprocessing --sample data/tsv/tiny-manta.tsv --tools Manta
nf_test . --step preprocessing --test --noReports

# Clean up images
if [ $PROFILE = travis ]
then
  docker rmi -f maxulysse/fastqc:1.1 maxulysse/mapreads:1.1 maxulysse/picard:1.1 maxulysse/runmanta:1.1
else
  rm -rf work/singularity/fastqc-1.1.img work/singularity/mapreads-1.1.img work/singularity/picard-1.1.img work/singularity/runmanta-1.1.img
fi

nf_test . --step skipPreprocessing --tools MuTect2,snpEff,VEP --noReports
nf_test . --step annotate --tools snpEff,VEP --annotateTools Strelka
nf_test . --step annotate --tools snpEff --annotateVCF VariantCalling/Manta/Manta_9876T_vs_1234N.diploidSV.vcf,VariantCalling/Manta/Manta_9876T_vs_1234N.somaticSV.vcf --noReports
