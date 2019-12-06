#!/bin/bash
set -xeuo pipefail

# This script run sarek tests
# https://github.com/nf-core/test-datasets/raw/sarek

usage() { echo "Usage: $0 <-p profile> <-t test> <-c cpus> <-n> <-v> <-m memory>" 1>&2; exit 1; }

CPUS=2
LOGS=''
MEMORY='6.GB'
NXF_SINGULARITY_CACHEDIR=${NXF_SINGULARITY_CACHEDIR:-work/singularity/.}
OFFLINE=false
PROFILE=docker
REPORTS=''
TEST=MULTIPLE
TRAVIS=${TRAVIS:-false}
TRAVIS_BUILD_DIR=${TRAVIS_BUILD_DIR:-.}
VERBOSE=''

while [[ $# -gt 0 ]]
do
  key=$1
  case $key in
    -c|--cpus)
    CPUS=$2
    shift # past value
    ;;
    -m|--memory)
    MEMORY=$2
    shift # past argument
    shift # past value
    ;;
    -n|--no-logs)
    LOGS=true
    shift # past value
    ;;
    --no-reports)
    REPORTS="--skipQC all"
    shift # past value
    ;;
    --offline)
    OFFLINE=true
    shift # past value
    ;;
    -p|--profile)
    PROFILE=$2
    shift # past argument
    shift # past value
    ;;
    -t|--test)
    TEST=$2
    shift # past argument
    shift # past value
    ;;
    -v|--verbose)
    VERBOSE="-ansi-log false -dump-channels"
    shift # past value
    ;;
    *) # unknown option
    usage
    shift # past argument
    ;;
  esac
done

function manage_logs() {
  if [[ $LOGS ]]
  then
    rm -rf .nextflow* results/ work/
  fi
}

function run_sarek() {
  nextflow run ${TRAVIS_BUILD_DIR}/main.nf -profile test,${PROFILE} ${VERBOSE} --monochrome_logs ${REPORTS} --max_memory ${MEMORY} $@
}

if [[ $OFFLINE == false ]]
then
  PATHTOSAMPLE="https://github.com/nf-core/test-datasets/raw/sarek/testdata"
  SUFFIX="-https"
else
  PATHTOSAMPLE="data/testdata"
  SUFFIX=""
fi

if [[ $TEST == "GERMLINE" ]] && [[ $OFFLINE == false ]]
then
  rm -rf data
  git clone --single-branch --branch sarek https://github.com/nf-core/test-datasets.git data
fi

case $TEST in
  ANNOTATEBOTH)
  ANNOTATOR="merge,snpEFF,VEP"
  TEST=ANNOTATE
  ;;
  ANNOTATESNPEFF)
  ANNOTATOR="snpEFF"
  TEST=ANNOTATE
  ;;
  ANNOTATEVEP)
  ANNOTATOR="VEP"
  TEST=ANNOTATE
  ;;
esac

case $TEST in
  GERMLINE)
  run_sarek --tools=false --input data/testdata/tiny/normal
  run_sarek --tools=false --input results/Preprocessing/TSV/duplicateMarked.tsv --step recalibrate  -resume
  run_sarek --tools HaplotypeCaller --input results/Preprocessing/TSV/recalibrated.tsv --step variantCalling
  ;;
  MINIMAL)
  run_sarek --tools manta,Mpileup,strelka --input ${PATHTOSAMPLE}/tsv/tiny-manta-normal${SUFFIX}.tsv --genome smallerGRCh37
  run_sarek --tools manta,Mpileup,strelka --input ${PATHTOSAMPLE}/tsv/tiny-manta-normal${SUFFIX}.tsv --genome smallerGRCh37 --noIntervals -resume
  run_sarek --tools manta,Mpileup,strelka --input ${PATHTOSAMPLE}/tsv/tiny-manta-normal${SUFFIX}.tsv --genome minimalGRCh37 -resume
  run_sarek --tools manta,Mpileup,strelka --input ${PATHTOSAMPLE}/tsv/tiny-manta-normal${SUFFIX}.tsv --genome minimalGRCh37 --noIntervals -resume
  ;;
  MULTIPLE)
  run_sarek --tools FreeBayes,HaplotypeCaller,Manta,Mpileup,Strelka,snpEff,VEP,merge --input ${PATHTOSAMPLE}/tsv/tiny-multiple${SUFFIX}.tsv
  ;;
esac

if [[ $TEST == "GERMLINE" ]] && [[ $OFFLINE == false ]]
then
  rm -rf data
fi
