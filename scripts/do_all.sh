#!/bin/bash
set -xeuo pipefail

PROFILE=singularity
PUSH=''
REPOSITORY=maxulysse
TAG=1.2.3
TOOL=docker

while [[ $# -gt 0 ]]
do
    key=$1
    case $key in
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

if [ $TOOL = docker ]
then
    nextflow run buildContainers.nf -profile ${PROFILE} --verbose --docker ${PUSH} --repository ${REPOSITORY} --tag ${TAG} --containers caw,fastqc,freebayes,gatk,igvtools,multiqc,mutect1,picard,qualimap,r-base,runallelecount,snpeff
    nextflow run buildContainers.nf -profile ${PROFILE} --verbose --docker ${PUSH} --repository ${REPOSITORY} --tag ${TAG} --containers snpeffgrch37,snpeffgrch38,vepgrch37,vepgrch38
else
    nextflow run buildContainers.nf -profile ${PROFILE} --verbose --singularity --repository ${REPOSITORY} --tag ${TAG} --containerPath containers/ --containers caw,fastqc,freebayes,gatk,igvtools,multiqc,mutect1,picard,qualimap,r-base,runallelecount,snpeffgrch37,snpeffgrch38,vepgrch37,vepgrch38
fi
