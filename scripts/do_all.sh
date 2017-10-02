#!/bin/bash
set -xeuo pipefail

PROFILE="singularity"
PUSH=""
REPOSITORY="--repository maxulysse"
TAG="1.1"
TOOL="docker"

while [[ $# -gt 1 ]]
do
    key="$1"
    case $key in
        -r|--repository)
        REPOSITORY="--repository $2"
        shift
        ;;
        -t|--tag)
        TAG="--tag $2"
        shift
        ;;
        -p|--profile)
        PROFILE="$2"
        shift
        ;;
        --push)
        PUSH=--push
        ;;
        --pull)
        TOOL=singularity
        ;;
        *) # unknown option
        ;;
    esac
    shift
done

if [ $TOOL = docker ]
then
    nextflow run buildContainers.nf -profile ${PROFILE} --verbose --docker ${PUSH} ${REPOSITORY} ${TAG} --containers caw,fastqc,freebayes,gatk,igvtools,multiqc,mutect1,picard,qualimap,runascat,runconvertallelecounts,snpeff,vep
    nextflow run buildContainers.nf -profile ${PROFILE} --verbose --docker ${PUSH} ${REPOSITORY} ${TAG} --containers runallelecount,snpeffgrch37,snpeffgrch38,vepgrch37,vepgrch38
else
    nextflow run buildContainers.nf -profile ${PROFILE} --verbose --singularity ${REPOSITORY} ${TAG} --singularityPublishDir containers/ --containers bcftools
    caw,fastqc,freebayes,gatk,igvtools,multiqc,mutect1,picard,qualimap,runallelecount,runascat,runconvertallelecounts,snpeffgrch37,snpeffgrch38,vepgrch37,vepgrch38
fi
