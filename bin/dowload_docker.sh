#!/bin/bash
set -xeuo pipefail

TEST=ALL
TRAVIS=${TRAVIS:-false}

TMPDIR=`pwd`/tmp
mkdir -p $TMPDIR

while [[ $# -gt 0 ]]
do
  key=$1
  case $key in
    -t|--test)
    TEST=$2
    shift # past argument
    shift # past value
    ;;
    *) # unknown option
    shift # past argument
    ;;
  esac
done

if [[ ALL,ANNOTATEALL,ANNOTATEVEP =~ $TEST ]]
then
  docker pull nfcore/sarekvep:dev.GRCh37
  docker tag nfcore/sarekvep:dev.GRCh37 nfcore/sarekvep:dev.GRCh37
elif [[ ALL,ANNOTATEALL,ANNOTATESNPEFF =~ $TEST ]]
then
  docker pull nfcore/sareksnpeff:dev.GRCh37
  docker tag nfcore/sareksnpeff:dev.GRCh37 nfcore/sareksnpeff:dev.GRCh37
else
  docker pull nfcore/sarek:dev
  docker tag nfcore/sarek:dev nfcore/sarek:dev
fi
