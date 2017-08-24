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

nf_test . --step preprocessing --test

# Clean up images
if [ $PROFILE = travis ]
then
  docker rmi -f maxulysse/fastqc:1.1 maxulysse/mapreads:1.1 maxulysse/picard:1.1
else
  rm -rf work/singularity/fastqc-1.1.img work/singularity/mapreads-1.1.img work/singularity/picard-1.1.img
fi

nf_test . --step realign --noReports
nf_test . --step realign --tools HaplotypeCaller
nf_test . --step realign --tools HaplotypeCaller --noReports --noGVCF
nf_test . --step recalibrate --noReports
nf_test . --step recalibrate --tools FreeBayes,HaplotypeCaller,MuTect1,MuTect2,Strelka

# Test whether restarting from an already recalibrated BAM works
nf_test . --step skipPreprocessing --tools Strelka --noReports
