#!/bin/bash
set -xeuo pipefail

ANNOTATE=false
ANNOTATEVCF=''
GENOME=GRCh38
GENOMEBASE=''
GERMLINE=false
PROFILE=singularity
SAMPLEDIR=''
SAMPLETSV=''
SOMATIC=false
STEP='mapping'
TOOLS='haplotypecaller,strelka,manta'
VARIANTCALLING=false

while [[ $# -gt 0 ]]
do
  key=$1
  case $key in
    -a|--annotate)
    ANNOTATE=true
    shift # past argument
    ;;
    -b|--genomeBase)
    GENOMEBASE=$2
    shift # past argument
    shift # past value
    ;;
    -c|--somatic)
    SOMATIC=true
    shift # past argument
    ;;
    -d|--sampleDir)
    SAMPLEDIR=$2
    shift # past argument
    shift # past value
    ;;
    -f|--annotateVCF)
    ANNOTATEVCF=$2
    shift # past argument
    shift # past value
    ;;
    -g|--genome)
    GENOME=$2
    shift # past argument
    shift # past value
    ;;
    -i|--sample)
    SAMPLETSV=$2
    shift # past argument
    shift # past value
    ;;
    -l|--germline)
    GERMLINE=true
    shift # past argument
    ;;
    -p|--profile)
    PROFILE=$2
    shift # past argument
    shift # past value
    ;;
    -s|--step)
    STEP=$2
    shift # past argument
    shift # past value
    ;;
    -t|--tools)
    TOOLS=$2
    shift # past argument
    shift # past value
    ;;
    -v|--variantCalling)
    TOOLS=$2
    shift # past argument
    shift # past value
    ;;
    *) # unknown option
    shift # past argument
    ;;
  esac
done

function run_sarek() {
  echo "$(tput setaf 1)nextflow run $@ -profile $PROFILE --genome $GENOME --genome_base $GENOMEBASE --verbose$(tput sgr0)"
  nextflow run $@ -profile $PROFILE --genome $GENOME --genome_base $GENOMEBASE --verbose
}

if [[ $GERMLINE == true ]] && [[ $SOMATIC == true ]]
then
  exit
fi

if [[ $GERMLINE == true ]] && [[ $ANNOTATE == true ]]
then
  exit
fi

if [[ $SAMPLEDIR != '' ]] && [[ $SAMPLETSV != '' ]]
then
  exit
fi

if [[ $SAMPLEDIR == '' ]] && [[ $SAMPLETSV == '' ]] && [[ $ANNOTATE == false ]]
then
  exit
fi

if [[ $SOMATIC == true ]] && [[ $SAMPLEDIR != '' ]]
then
  exit
fi

if [[ $GERMLINE == true ]] && [[ $SAMPLEDIR != '' ]]
then
  run_sarek main.nf --step $STEP --sampleDir $SAMPLEDIR
fi

if [[ $GERMLINE == true ]] && [[ $SAMPLETSV != '' ]]
then
  run_sarek main.nf --step $STEP --sample $SAMPLETSV
fi

if [[ $GERMLINE == true ]] && [[ $VARIANTCALLING == true ]]
then
  run_sarek germlineVC.nf --tools $TOOLS
fi

if [[ $SOMATIC == true ]] && [[ $SAMPLETSV != '' ]]
then
  run_sarek main.nf --step $STEP --sample $SAMPLETSV
fi

if [[ $SOMATIC == true ]] && [[ $VARIANTCALLING == true ]]
then
  run_sarek germlineVC.nf --tools $TOOLS
  run_sarek somaticVC.nf --tools $TOOLS
fi


if [[ $ANNOTATE == true ]]
then
  run_sarek annotate.nf --tools $TOOLS --annotateVCF $ANNOTATEVCF
fi
