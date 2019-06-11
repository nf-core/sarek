#!/bin/bash
set -euo pipefail

TEST=false
NAME=sarek-$(git describe --tags --always)

while [[ $# -gt 0 ]]
do
  key="$1"
  case $key in
    -i|-t|--include-test-data)
    TEST=true
    shift # past argument
    ;;
    *) # unknown option
    shift # past argument
    ;;
  esac
done

echo "Archiving nf-core/configs"
rm -rf configs
git submodule add -f https://github.com/nf-core/configs.git configs

if [[ $TEST == false ]]
then
  echo "Archiving nf-core/sarek"
  git-archive-all --prefix=${NAME} --force-submodules ${NAME}.tar.gz
else
  echo "Archiving nf-core/test-datasets:sarek"
  rm -rf data
  git submodule add -f --branch sarek https://github.com/nf-core/test-datasets.git data
  echo "Archiving nf-core/sarek"
  git-archive-all --prefix=${NAME} --force-submodules ${NAME}.tar.gz
fi

echo "Wrote ${NAME}.tar.gz"
