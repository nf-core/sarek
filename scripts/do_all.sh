#!/bin/bash
set -xeuo pipefail

PROFILE="singularity"
PUSH=""
REPOSITORY="--repository maxulysse"
TAG="1.2.1"
TOOL="docker"

while [[ $# -gt 0 ]]
do
    key="$1"
    case $key in
        -r|--repository)
        REPOSITORY="--repository $2"
        shift # past argument
        shift # past value
        ;;
        -t|--tag)
        TAG="--tag $2"
        shift # past argument
        shift # past value
        ;;
        -p|--profile)
        PROFILE="$2"
        shift # past argument
        shift # past value
        ;;
        --push)
        PUSH=--push
        shift # past argument
        ;;
        --pull)
        TOOL=singularity
        shift # past argument
        ;;
        *) # unknown option
        shift # past argument
        ;;
    esac
done

if [ $TOOL = docker ]
then
    nextflow run buildContainers.nf -profile ${PROFILE} --verbose --docker ${PUSH} ${REPOSITORY} ${TAG} --containers caw,fastqc,freebayes,gatk,igvtools,multiqc,mutect1,picard,qualimap,r-base,runallelecount,snpeff
    nextflow run buildContainers.nf -profile ${PROFILE} --verbose --docker ${PUSH} ${REPOSITORY} ${TAG} --containers snpeffgrch37,snpeffgrch38,vepgrch37,vepgrch38
else
    nextflow run buildContainers.nf -profile ${PROFILE} --verbose --singularity ${REPOSITORY} ${TAG} --singularityPublishDir containers/ --containers caw,fastqc,freebayes,gatk,igvtools,multiqc,mutect1,picard,qualimap,r-base,runallelecount,snpeffgrch37,snpeffgrch38,vepgrch37,vepgrch38
fi
