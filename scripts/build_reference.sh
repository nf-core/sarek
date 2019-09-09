#!/bin/bash
set -xeuo pipefail

# This script build small reference for sarek tests
# https://github.com/nf-core/test-datasets/raw/sarek

usage() { echo "Usage: $0 <-p profile> <-t test> <-v> <-m memory>" 1>&2; exit 1; }

MEMORY='6.GB'
NXF_SINGULARITY_CACHEDIR=${NXF_SINGULARITY_CACHEDIR:-work/singularity/.}
OFFLINE=''
PROFILE=docker
TEST=ALL
TRAVIS=${TRAVIS:-false}
TRAVIS_BUILD_DIR=${TRAVIS_BUILD_DIR:-.}
VERBOSE=''

while [[ $# -gt 0 ]]
do
  key=$1
  case $key in
    -t|--test)
    TEST=$2
    shift # past argument
    shift # past value
    ;;
    -m|--memory)
    MEMORY=$2
    shift # past argument
    shift # past value
    ;;
    --offline)
    OFFLINE="--offline"
    shift # past value
    ;;
    -p|--profile)
    PROFILE=$2
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

# Build references for smallGRCh37
if ! [[ ANNOTATEBOTH,ANNOTATESNPEFF,ANNOTATEVEP,LINT =~ $TEST ]]
 then
  rm -rf references
  nextflow run ${TRAVIS_BUILD_DIR}/build.nf -profile test,${PROFILE} --build --outdir references ${VERBOSE} ${OFFLINE} --max_memory ${MEMORY}
  rm -rf .nextflow* references/pipeline_info work
fi
