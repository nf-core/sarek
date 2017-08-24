#!/bin/bash
set -xeuo pipefail
PROFILE=$1
TEST=$2

function nf_test() {
  echo "$(tput setaf 1)nextflow run $@ -profile $PROFILE -resume --verbose$(tput sgr0)"
  nextflow run "$@" -profile $PROFILE -resume --verbose
}

nf_test buildReferences.nf --download

if [ $TEST = MAPPING ]
then
  nf_test . --test --step preprocessing
fi

if [ $TEST = REALIGN ]
then
  nf_test . --test --step preprocessing
  nf_test . --step realign --noReports
  nf_test . --step realign --tools HaplotypeCaller
  nf_test . --step realign --tools HaplotypeCaller --noReports --noGVCF
fi

if [ $TEST = RECALIBRATE ]
then
  nf_test . --test --step preprocessing
  nf_test . --step recalibrate --noReports
  nf_test . --step recalibrate --tools FreeBayes,HaplotypeCaller,MuTect1,MuTect2,Strelka
  # Test whether restarting from an already recalibrated BAM works
  nf_test . --step skipPreprocessing --tools Strelka --noReports
fi

if [ $TEST = ANNOTATE ]
then
  nf_test . --step preprocessing --sample data/tsv/tiny-manta.tsv --tools Manta
  nf_test . --test --step preprocessing --tools MuTect2,Strelka,snpEff,VEP --noReports
  nf_test . --step annotate --tools snpEff,VEP --annotateTools Strelka
  nf_test . --step annotate --tools snpEff --annotateVCF VariantCalling/Manta/Manta_9876T_vs_1234N.diploidSV.vcf,VariantCalling/Manta/Manta_9876T_vs_1234N.somaticSV.vcf --noReports
fi
