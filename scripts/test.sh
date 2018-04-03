#!/bin/bash
set -xeuo pipefail

BUILD=false
GENOME=smallGRCh37
PROFILE=singularity
SAMPLE=data/tsv/tiny.tsv
TEST=ALL
TRAVIS=${TRAVIS:-false}

while [[ $# -gt 0 ]]
do
  key=$1
  case $key in
    -g|--genome)
    GENOME=$2
    shift # past argument
    shift # past value
    ;;
    -p|--profile)
    PROFILE=$2
    shift # past argument
    shift # past value
    ;;
    -s|--sample)
    SAMPLE=$2
    shift # past argument
    shift # past value
    ;;
    -t|--test)
    TEST=$2
    shift # past argument
    shift # past value
    ;;
    -b|--build)
    BUILD=true
    shift # past value
    ;;
    *) # unknown option
    shift # past argument
    ;;
  esac
done

function nf_test() {
  echo "$(tput setaf 1)nextflow run $@ -profile $PROFILE --genome $GENOME -resume --verbose$(tput sgr0)"
  nextflow run $@ -profile $PROFILE --genome $GENOME -resume --genome_base $PWD/References/$GENOME --verbose
}

function run_wrapper() {
  ./scripts/wrapper.sh $@ --profile $PROFILE --genome $GENOME --genomeBase $PWD/References/$GENOME --verbose
}

function clean_repo() {
  if [[ $TRAVIS == false ]]
  then
    rm -rf work .nextflow* Preprocessing Reports Annotation VariantCalling Results
  fi
}

# Build references only for smallGRCh37
if [[ $GENOME == smallGRCh37 ]] && [[ $TEST != BUILDCONTAINERS ]] && [[ BUILD ]]
then
  nf_test buildReferences.nf --download --outDir References/$GENOME
  # Remove images only on TRAVIS
  if [[ $PROFILE == docker ]] && [[ $TRAVIS == true ]]
  then
    docker rmi -f maxulysse/igvtools:latest
  elif [[ $PROFILE == singularity ]] && [[ $TRAVIS == true ]]
  then
    rm -rf work/singularity/igvtools-latest.img
  fi
fi

if [[ ALL,DIR =~ $TEST ]]
then
  run_wrapper --germline --sampleDir data/tiny/tiny/normal
  clean_repo
fi

if [[ ALL,STEP =~ $TEST ]]
then
  run_wrapper --germline --sample $SAMPLE
  run_wrapper --germline --step realign --noReports
  run_wrapper --germline --step recalibrate --noReports
  clean_repo
fi

if [[ ALL,GERMLINE =~ $TEST ]]
then
  run_wrapper --germline --sampleDir data/tiny/tiny/normal --variantCalling --tools HaplotypeCaller
  clean_repo
fi

if [[ ALL,TOOLS =~ $TEST ]]
then
  run_wrapper --somatic --sample $SAMPLE --variantCalling  --tools FreeBayes,HaplotypeCaller,MuTect1,MuTect2
fi

if [[ ALL,MANTA =~ $TEST ]]
then
  run_wrapper --somatic --sample data/tsv/tiny-manta.tsv --variantCalling --tools Manta --noReports
  run_wrapper --somatic --sample data/tsv/tiny-manta.tsv --variantCalling --tools Manta,Strelka --noReports --strelkaBP
  clean_repo
fi


if [[ ALL,ANNOTATESNPEFF,ANNOTATEVEP =~ $TEST ]]
then
  if [[ $TEST = ANNOTATESNPEFF ]]
  then
    ANNOTATOR=snpEFF
  elif [[ $TEST = ANNOTATEVEP ]]
  then
    ANNOTATOR=VEP
  elif  [[ $TEST = ALL ]]
  then
    ANNOTATOR=snpEFF,VEP
  fi
  if [[ $PROFILE == docker ]] && [[ $TRAVIS == true ]]
  then
    docker rmi -f maxulysse/sarek:latest
    docker rmi -f maxulysse/picard:latest
  elif [[ $PROFILE == singularity ]] && [[ $TRAVIS == true ]]
  then
    rm -rf work/singularity/sarek-latest.img
    rm -rf work/singularity/picard-latest.img
  fi
  run_wrapper --annotate --tools ${ANNOTATOR} --annotateVCF data/tiny/vcf/Strelka_1234N_variants.vcf.gz --noReports
  run_wrapper --annotate --tools ${ANNOTATOR} --annotateVCF data/tiny/vcf/Strelka_1234N_variants.vcf.gz,data/tiny/vcf/Strelka_9876T_variants.vcf.gz --noReports
  clean_repo
fi

if [[ ALL,BUILDCONTAINERS =~ $TEST ]] && [[ $PROFILE == docker ]]
then
  nf_test buildContainers.nf --docker --containers fastqc,gatk,igvtools,multiqc,mutect1,picard,qualimap,runallelecount,r-base,snpeff,sarek
  clean_repo
fi
