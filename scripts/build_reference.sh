#!/bin/bash
set -xeuo pipefail

# This script build small reference for sarek tests
# https://github.com/nf-core/test-datasets/raw/sarek

usage() { echo "Usage: $0 <-p profile> <-t test> <-v>" 1>&2; exit 1; }

NXF_SINGULARITY_CACHEDIR=${NXF_SINGULARITY_CACHEDIR:-work/singularity/.}
OFFLINE=''
PROFILE=docker
TEST=ALL
TRAVIS_BUILD_DIR=${TRAVIS_BUILD_DIR:-.}
TRAVIS=${TRAVIS:-false}
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
if [[ $TEST != ANNOTATESNPEFF ]] && [[ $TEST != ANNOTATEVEP ]]
then
  rm -rf references
  nextflow run ${TRAVIS_BUILD_DIR}/build.nf -profile test,${PROFILE} --build --outdir references ${VERBOSE} ${OFFLINE}
  rm -rf .nextflow* references/pipeline_info work
fi
