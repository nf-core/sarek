#!/bin/bash
set -xeuo pipefail

CPUS=2
TEST=ALL
TRAVIS_BUILD_DIR=${TRAVIS_BUILD_DIR:-.}
TRAVIS=${TRAVIS:-false}

while [[ $# -gt 0 ]]
do
  key=$1
  case $key in
    -t|--test)
    TEST=$2
    shift # past argument
    shift # past value
    ;;
    -c|--cpus)
    CPUS=$2
    shift # past value
    ;;
    *) # unknown option
    shift # past argument
    ;;
  esac
done

function run_sarek() {
  nextflow run ${TRAVIS_BUILD_DIR}/main.nf -profile test,docker -ansi-log false -dump-channels $@
}

if [[ ALL,GERMLINE =~ $TEST ]]
then
  rm -rf data
  git clone --single-branch --branch sarek https://github.com/nf-core/test-datasets.git data
  run_sarek --sample data/testdata/tiny/normal --tools HaplotypeCaller,Strelka --noReports
  run_sarek --step recalibrate --sample results/Preprocessing/TSV/duplicateMarked.tsv --noReports
fi

if [[ ALL,SOMATIC =~ $TEST ]]
then
	run_sarek --tools FreeBayes,HaplotypeCaller,Manta,Strelka,Mutect2 --noReports
fi

if [[ ALL,TARGETED =~ $TEST ]]
then
	run_sarek --tools FreeBayes,HaplotypeCaller,Manta,Strelka,Mutect2 --noReports --targetBED https://github.com/nf-core/test-datasets/raw/sarek/testdata/target.bed
fi

if [[ ALL,ANNOTATEALL,ANNOTATESNPEFF,ANNOTATEVEP =~ $TEST ]]
then
  if [[ $TEST = ANNOTATESNPEFF ]]
  then
    ANNOTATOR=snpEFF
  elif [[ $TEST = ANNOTATEVEP ]]
  then
    ANNOTATOR=VEP
  elif [[ ALL,ANNOTATEALL =~ $TEST ]]
  then
    ANNOTATOR=merge,snpEFF,VEP
  fi
  run_sarek --step annotate --tools ${ANNOTATOR} --sample https://github.com/nf-core/test-datasets/raw/sarek/testdata/vcf/Strelka_1234N_variants.vcf.gz --noReports
fi

if [[ MULTIPLE =~ $TEST ]]
then
  run_sarek --sample https://github.com/nf-core/test-datasets/raw/sarek/testdata/tsv/tiny-multiple.tsv --tools FreeBayes,HaplotypeCaller,Manta,Strelka,Mutect2 --noReports
fi
