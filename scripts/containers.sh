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
  elif [[ $TEST = ANNOTATESNPEFF ]]
  then
    docker pull maxulysse/snpeffgrch37:latest
  fi
  docker pull maxulysse/sarek:latest
fi
