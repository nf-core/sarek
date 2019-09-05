#!/bin/bash
set -xeuo pipefail

# This script download and tag image for sarek tests

usage() { echo "Usage: $0 <-t test|annotation tool> <-n engine> <-T version to pull/build> <-g genome>" 1>&2; exit 1; }

ENGINE=docker
GENOME=smallGRCh37
NXF_SINGULARITY_CACHEDIR=${NXF_SINGULARITY_CACHEDIR:-work/singularity/.}
TEST=ALL
VERSION=dev

while [[ $# -gt 0 ]]
do
  key=$1
  case $key in
    -g|--genome)
    GENOME=$2
    shift # past argument
    shift # past value
    ;; 
    -n|--engine)
    ENGINE=$2
    shift # past argument
    shift # past value
    ;;
    -T|--tagged-version)
    VERSION=$2
    shift # past argument
    shift # past value
    ;;
    -t|--test|--tool) 
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

SOURCEGENOME=${GENOME}

if [[ smallGRCh37 =~ $SOURCEGENOME ]]
then
  SOURCEGENOME=GRCh37
fi

get_image(){
  CONTAINER=$1
  SOURCE=$2
  TARGET=$3
  if [[ docker =~ $ENGINE ]]
  then
    docker pull nfcore/${1}:${2}
    docker tag nfcore/${1}:${2} nfcore/${1}:${3}
  elif  [[ singularity =~ $ENGINE ]]
  then
    mkdir -p ${NXF_SINGULARITY_CACHEDIR}
    singularity build ${NXF_SINGULARITY_CACHEDIR}/nfcore-${1}-${3}.img docker://nfcore/${1}:${2}
  fi
}

if [[ ALL,ANNOTATEBOTH,ANNOTATESNPEFF,SNPEFF =~ $TEST ]]
then
  get_image sareksnpeff ${VERSION}.${SOURCEGENOME} ${VERSION}.${GENOME}
fi

if [[ ALL,ANNOTATEBOTH,ANNOTATEVEP,VEP =~ $TEST ]]
then
  get_image sarekvep ${VERSION}.${SOURCEGENOME} ${VERSION}.${GENOME}
fi

if ! [[ ANNOTATEBOTH,ANNOTATESNPEFF,ANNOTATEVEP,LINT,SNPEFF,VEP =~ $TEST ]]
then
  get_image sarek ${VERSION} ${VERSION}
fi
