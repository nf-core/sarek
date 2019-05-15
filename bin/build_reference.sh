#!/bin/bash
set -xeuo pipefail

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
    shift # past argument
    ;;
  esac
done

# Build references for smallGRCh37
if [[ $TEST != ANNOTATESNPEFF ]] && [[ $TEST != ANNOTATEVEP ]]
then
  rm -rf references
  nextflow run ${TRAVIS_BUILD_DIR}/build.nf -profile test,${PROFILE} --build --outdir references ${VERBOSE}
  rm -rf .nextflow* references/pipeline_info work
fi
