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
  nextflow run ${TRAVIS_BUILD_DIR}/main.nf -profile docker -ansi-log false --publishDirMode link --max_memory 7.GB --max_cpus 2 -dump-channels --genome smallGRCh37 --igenomes_base reference $@
}

if [[ ALL,GERMLINE =~ $TEST ]]
then
  run_sarek --sample data/testdata/tiny/normal --tools HaplotypeCaller,Strelka --noReports
  run_sarek --step recalibrate --noReports
	clean_repo
fi

if [[ ALL,SOMATIC =~ $TEST ]]
then
	run_sarek --sample data/testdata/tsv/tiny-manta.tsv --tools FreeBayes,HaplotypeCaller,Manta,Strelka,Mutect2 --noReports
	clean_repo
fi

if [[ ALL,TARGETED =~ $TEST ]]
then
	run_sarek --sample data/testdata/tsv/tiny-manta.tsv --tools FreeBayes,HaplotypeCaller,Manta,Strelka,Mutect2 --noReports --targetBED data/testdata/target.bed
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
  run_sarek --step annotate --tools ${ANNOTATOR} --annotateVCF data/testdata/vcf/Strelka_1234N_variants.vcf.gz --noReports
  clean_repo
fi

if [[ MULTIPLE =~ $TEST ]]
then
  run_sarek --sample data/testdata/tsv/tiny-multiple.tsv --tools FreeBayes,HaplotypeCaller,Manta,Strelka,Mutect2 --noReports
fi
