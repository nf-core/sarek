#!/bin/bash
set -euo pipefail

# This script makes an archive of sarek, with or without configs and test datasets
# https://github.com/nf-core/sarek

usage() { echo "Usage: $0 <-t> <-c>" 1>&2; exit 1; }

CONFIGS=false
NAME=sarek-$(git describe --tags --always)
TEST=false

while [[ $# -gt 0 ]]
do
  key="$1"
  case $key in
    -i|-t|--include-test-data)
    TEST=true
    shift # past argument
    ;;
    -c|--include-configs)
    CONFIGS=true
    shift # past argument
    ;;
    *) # unknown option
    shift # past argument
    usage
    ;;
  esac
done

if [[ $CONFIGS == true ]]
then
  echo "Archiving nf-core/configs"
  git submodule add -f https://github.com/nf-core/configs.git configs
fi

if [[ $TEST == true ]]
then
  echo "Archiving nf-core/test-datasets:sarek"
  git submodule add -f --branch sarek https://github.com/nf-core/test-datasets.git data
fi

echo "Archiving nf-core/sarek"

if [[ $CONFIGS == true ]] || [[ $TEST == true ]]
then
  git-archive-all --prefix=${NAME} --force-submodules ${NAME}.tar.gz
else
  git archive --format=tar.gz HEAD --prefix=${NAME}/ > ${NAME}.tar.gz
fi

echo "Wrote ${NAME}.tar.gz"
