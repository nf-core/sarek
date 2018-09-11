#!/bin/bash
set -xeuo pipefail

PROFILE=singularity
TEST=ALL
TRAVIS=${TRAVIS:-false}

TMPDIR=`pwd`/tmp
mkdir -p $TMPDIR

while [[ $# -gt 0 ]]
do
  key=$1
  case $key in
    -p|--profile)
    PROFILE=$2
    shift # past argument
    shift # past value
    ;;
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

if [[ $PROFILE = docker ]] && [[ $TRAVIS == true ]]
then
  if [[ $TEST = ANNOTATEVEP ]]
  then
    docker pull maxulysse/vepgrch37:latest
  else
    docker pull maxulysse/snpeffgrch37:latest
  fi
fi

if [[ $TEST = ANNOTATESNPEFF ]] && [[ $PROFILE = singularity ]] && [[ $TRAVIS == true ]]
then
  singularity build $TMPDIR/maxulysse-snpeffgrch37-latest.simg docker://maxulysse/snpeffgrch37:latest
fi
