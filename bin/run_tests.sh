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

if [[ ALL,GERMLINE =~ $TEST ]]
then
  nextflow run ${TRAVIS_BUILD_DIR}/main.nf --sample data/testdata/tiny/normal --tools HaplotypeCaller,Strelka --noReports
  nextflow run ${TRAVIS_BUILD_DIR}/main.nf --step recalibrate --noReports
	clean_repo
fi

if [[ ALL,SOMATIC =~ $TEST ]]
then
	nextflow run ${TRAVIS_BUILD_DIR}/main.nf --sample data/testdata/tsv/tiny-manta.tsv --tools FreeBayes,HaplotypeCaller,Manta,Strelka,Mutect2 --noReports
	clean_repo
fi

if [[ ALL,TARGETED =~ $TEST ]]
then
	nextflow run ${TRAVIS_BUILD_DIR}/main.nf --sample data/testdata/tsv/tiny-manta.tsv --tools FreeBayes,HaplotypeCaller,Manta,Strelka,Mutect2 --noReports --targetBED data/testdata/target.bed
	clean_repo
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
  nextflow run ${TRAVIS_BUILD_DIR}/main.nf --step annotate --tools ${ANNOTATOR} --annotateVCF data/testdata/vcf/Strelka_1234N_variants.vcf.gz --noReports
  clean_repo
fi

if [[ MULTIPLE =~ $TEST ]]
then
  nextflow run ${TRAVIS_BUILD_DIR}/main.nf --sample data/testdata/tsv/tiny-multiple.tsv --tools FreeBayes,HaplotypeCaller,Manta,Strelka,Mutect2 --noReports
fi
