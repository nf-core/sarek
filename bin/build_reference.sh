#!/bin/bash
set -xeuo pipefail

BUILD=false
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
    -b|--build)
    BUILD=true
    shift # past value
    ;;
    *) # unknown option
    shift # past argument
    ;;
  esac
done

# Always download test data
rm -rf data
git clone --single-branch --branch sarek https://github.com/nf-core/test-datasets.git data

# Build references for smallGRCh37
if [[ BUILD ]] && [[ $TEST != ANNOTATESNPEFF ]] && [[ $TEST != ANNOTATEVEP ]]
then
  rm -rf references
  nextflow run build.nf -profile docker --genome smallGRCh37 --refdir data/reference --outdir references --publishDirMode link --max_memory 7.GB --max_cpus 2 -ansi-log false
  rm -rf .nextflow* references/pipeline_info work
fi
