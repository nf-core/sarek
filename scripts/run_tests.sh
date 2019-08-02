#!/bin/bash
set -xeuo pipefail

# This script run sarek tests
# https://github.com/nf-core/test-datasets/raw/sarek

usage() { echo "Usage: $0 <-p profile> <-t test> <-c cpus> <-n> <-v>" 1>&2; exit 1; }

CPUS=2
LOGS=''
REPORTS=''
NXF_SINGULARITY_CACHEDIR=${NXF_SINGULARITY_CACHEDIR:-work/singularity/.}
OFFLINE=false
PROFILE=docker
TEST=ALL
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
    -n|--no-logs)
    LOGS=true
    shift # past value
    ;;
    --no-reports)
    REPORTS="--skip all"
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
  nextflow run ${TRAVIS_BUILD_DIR}/main.nf -profile test,${PROFILE} ${VERBOSE} --monochrome_logs ${REPORTS} $@
}

if [[ ALL,GERMLINE =~ $TEST ]]
then
  if [[ $OFFLINE == false ]]
  then
    rm -rf data
    git clone --single-branch --branch sarek https://github.com/nf-core/test-datasets.git data
  fi
  run_sarek --tools=false --sample data/testdata/tiny/normal
  run_sarek --tools=false --sample results/Preprocessing/TSV/duplicateMarked.tsv --step recalibrate
  run_sarek --tools HaplotypeCaller --sample results/Preprocessing/TSV/recalibrated.tsv --step variantCalling
  if [[ $OFFLINE == false ]]
  then
    rm -rf data
  fi
  manage_logs
fi

if [[ ALL,SOMATIC =~ $TEST ]]
then
  OPTIONS="--tools FreeBayes,HaplotypeCaller,Manta,Strelka,TIDDIT,Mutect2"
  if [[ $OFFLINE == false ]]
  then
    run_sarek ${OPTIONS}
  else
    run_sarek ${OPTIONS} --sample data/testdata/tsv/tiny-manta.tsv
  fi
  manage_logs
fi

if [[ ALL,TARGETED =~ $TEST ]]
then
  OPTIONS="--tools FreeBayes,HaplotypeCaller,Manta,Strelka,TIDDIT,Mutect2"
  if [[ $OFFLINE == false ]]
  then
    run_sarek ${OPTIONS} --targetBED https://github.com/nf-core/test-datasets/raw/sarek/testdata/target.bed
  else
    run_sarek ${OPTIONS} --sample data/testdata/tsv/tiny-manta.tsv --targetBED data/testdata/target.bed
  fi
  manage_logs
fi

if [[ $OFFLINE == false ]]
then
  pathToSample="https://github.com/nf-core/test-datasets/raw/sarek/testdata"
else
  pathToSample="data/testdata"
fi

if [[ ALL,ANNOTATEBOTH,ANNOTATESNPEFF,ANNOTATEVEP =~ $TEST ]]
then
  if [[ $TEST = ANNOTATESNPEFF ]]
  then
    ANNOTATOR=snpEFF
  elif [[ $TEST = ANNOTATEVEP ]]
  then
    ANNOTATOR=VEP
  elif [[ ALL,ANNOTATEBOTH =~ $TEST ]]
  then
    ANNOTATOR=merge,snpEFF,VEP
  fi
  run_sarek --step annotate --tools ${ANNOTATOR} --sample ${pathToSample}/vcf/Strelka_1234N_variants.vcf.gz
  manage_logs
fi

if [[ MULTIPLE =~ $TEST ]]
then
  OPTIONS="--tools FreeBayes,HaplotypeCaller,Manta,Strelka,TIDDIT,Mutect2,snpEff,VEP,merge"
  if [[ $OFFLINE == false ]]
  then
    run_sarek ${OPTIONS} --sample https://github.com/nf-core/test-datasets/raw/sarek/testdata/tsv/tiny-multiple-https.tsv
  else
    run_sarek ${OPTIONS} --sample data/testdata/tsv/tiny-multiple.tsv
  fi
  manage_logs
fi
