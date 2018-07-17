#!/bin/bash
set -euo pipefail

TEST=false
NAME=Sarek-$(git describe --tags)

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

if [[ $TEST == false ]]
then
  echo "Archiving Sarek without test data"
  git archive HEAD --prefix=$NAME/ | gzip > $NAME.tar.gz
else
  echo "Archiving Sarek with test data"
  git-archive-all --prefix=$NAME/ --force-submodules $NAME.tar.gz
fi

echo "Wrote $NAME.tar.gz"
