#!/bin/bash
set -xeuo pipefail

PROFILE=singularity
PUSH=''
REPOSITORY=maxulysse
GENOME=GRCh38
TAG=latest
TOOL=docker

while [[ $# -gt 0 ]]
do
    key=$1
    case $key in
        --genome)
        GENOME=$2
        shift # past argument
        shift # past value
        ;;
        -p|--profile)
        PROFILE=$2
        shift # past argument
        shift # past value
        ;;
        --pull)
        TOOL=singularity
        shift # past argument
        ;;
        --push)
        PUSH=--push
        shift # past argument
        ;;
        -r|--repository)
        REPOSITORY=$2
        shift # past argument
        shift # past value
        ;;
        -t|--tag)
        TAG=$2
        shift # past argument
        shift # past value
        ;;
        *) # unknown option
        shift # past argument
        ;;
    esac
done

if [[ $GENOME = smallGRCh37 ]]
then
    GENOME=GRCh37
fi

function toLower() {
    echo $1 | tr '[:upper:]' '[:lower:]'
}

if [[ $TOOL = docker ]] && [[ GRCh37,GRCh38 =~ $GENOME ]]
then
    nextflow run build.nf -profile ${PROFILE} -dump-channels --docker ${PUSH} --repository ${REPOSITORY} --tag ${TAG} --containers r-base,runallelecount,sarek,snpeff$(toLower ${GENOME}),vep$(toLower ${GENOME})
else
    nextflow run build.nf -profile ${PROFILE} -dump-channels --singularity --repository ${REPOSITORY} --tag ${TAG} --containerPath containers/ --containers r-base,runallelecount,sarek,snpeff$(toLower ${GENOME}),vep$(toLower ${GENOME})
fi
