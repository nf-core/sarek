#!/bin/bash
set -xeuo pipefail

PROFILE=singularity
TEST=ALL
TRAVIS=${TRAVIS:-false}

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

if [[ $TEST = ANNOTATEVEP ]] && [[ $PROFILE = docker ]] && [[ $TRAVIS == true ]]
then
  docker pull maxulysse/vepgrch37:latest
fi
