#!/bin/bash
set -xeuo pipefail

PUSH=""
REPOSITORY="maxulysse"
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
    nextflow run BuildContainers.nf --docker ${PUSH} ${REPOSITORY} ${TAG} --containers
    caw,fastqc,gatk,igvtools,multiqc,mutect1,picard,qualimap,runascat,runconvertallelecounts,snpeff,vep
    nextflow run BuildContainers.nf --docker ${PUSH} ${REPOSITORY} ${TAG} --containers runallelecount,snpeffgrch37,snpeffgrch38,vepgrch37,vepgrch38
else
    nextflow run BuildContainers.nf --singularity ${REPOSITORY} ${TAG} --singularityPublishDir containers/ --containers bcftools
    caw,fastqc,gatk,igvtools,multiqc,mutect1,picard,qualimap,runallelecount,runascat,runconvertallelecounts,snpeffgrch37,snpeffgrch38,vepgrch37,vepgrch38
fi
