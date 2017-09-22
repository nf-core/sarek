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
    nextflow run . --docker ${PUSH} ${REPOSITORY} ${TAG} --containers bcftools,concatvcf,fastqc,gatk,htslib,igvtools,multiqc,mutect1,picard,qualimap,runascat,runconvertallelecounts,samtools,strelka,snpeff,vep
    nextflow run . --docker ${PUSH} ${REPOSITORY} ${TAG} --containers mapreads,runallelecount,runmanta,snpeffgrch37,snpeffgrch38,vepgrch37,vepgrch38
else
    nextflow run . --singularity ${REPOSITORY} ${TAG} --singularityPublishDir containers/ --containers bcftools,concatvcf,fastqc,gatk,htslib,igvtools,mapreads,multiqc,mutect1,picard,qualimap,runallelecount,runascat,runconvertallelecounts,runmanta,samtools,snpeffgrch37,snpeffgrch38,strelka,vepgrch37,vepgrch38
fi
