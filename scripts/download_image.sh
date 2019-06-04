#!/bin/bash
set -xeuo pipefail

# This script download and tag image for sarek tests

usage() { echo "Usage: $0 <-t test> <-n engine>" 1>&2; exit 1; }

ENGINE=docker
NXF_SINGULARITY_CACHEDIR=${NXF_SINGULARITY_CACHEDIR:-work/singularity/.}
TEST=ALL

while [[ $# -gt 0 ]]
do
  key=$1
  case $key in
    -n|--engine)
    ENGINE=$2
    shift # past argument
    shift # past value
    ;;
    -t|--test)
    TEST=$2
    shift # past argument
    shift # past value
    ;;
    *) # unknown option
    usage
    shift # past argument
    ;;
  esac
done

if [[ ALL,ANNOTATEALL,ANNOTATESNPEFF =~ $TEST ]]
then
  if [[ docker =~ $ENGINE ]]
  then
    docker pull nfcore/sareksnpeff:dev.GRCh37
    docker tag nfcore/sareksnpeff:dev.GRCh37 nfcore/sareksnpeff:dev.smallGRCh37
  elif [[ singularity =~ $ENGINE ]]
  then
    mkdir -p work/singularity
    singularity build nfcore-sareksnpeff-dev.GRCh37.img docker://nfcore/sareksnpeff:dev.GRCh37
    mv nfcore-sareksnpeff-dev.GRCh37.img ${NXF_SINGULARITY_CACHEDIR}/.
  fi
fi

if [[ ALL,ANNOTATEALL,ANNOTATEVEP =~ $TEST ]]
then
  if [[ docker =~ $ENGINE ]]
  then
    docker pull nfcore/sarekvep:dev.GRCh37
    docker tag nfcore/sarekvep:dev.GRCh37 nfcore/sarekvep:dev.smallGRCh37
  elif [[ singularity =~ $ENGINE ]]
  then
    mkdir -p work/singularity
    singularity build nfcore-sarekvep-dev.GRCh37.img docker://nfcore/sarekvep:dev.GRCh37
    mv nfcore-sarekvep-dev.GRCh37.img ${NXF_SINGULARITY_CACHEDIR}/.
  fi
fi

if [[ ANNOTATEALL,ANNOTATEVEP,ANNOTATESNPEFF != $TEST ]]
then
  if [[ docker =~ $ENGINE ]]
  then
    docker pull nfcore/sarek:dev
    docker tag nfcore/sarek:dev nfcore/sarek:dev
  elif [[ singularity =~ $ENGINE ]]
  then
    mkdir -p work/singularity
    singularity build nfcore-sarek-dev.img docker://nfcore/sarek:dev
    mv nfcore-sarek-dev.img ${NXF_SINGULARITY_CACHEDIR}/.
  fi
fi
