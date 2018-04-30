#!/bin/bash
set -xeuo pipefail

BUILD=false
KEEP=false
GENOME=smallGRCh37
PROFILE=singularity
SAMPLE=Sarek-data/testdata/tsv/tiny.tsv
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
    -k|--keep)
    KEEP=true
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

function run_wrapper() {
  ./scripts/wrapper.sh $@ --profile $PROFILE --genome $GENOME --genomeBase $PWD/References/$GENOME --verbose
}

function clean_repo() {
  if [[ $TRAVIS == false ]] && [[ $KEEP == false ]]
  then
    echo "$(tput setaf 1)Cleaning directory$(tput sgr0)"
    rm -rf work .nextflow* Preprocessing Reports Annotation VariantCalling Results
  fi
}

# Build references only for smallGRCh37
if [[ $GENOME == smallGRCh37 ]] && [[ $TEST != BUILDCONTAINERS ]] && [[ BUILD ]]
then
  if [[ ! -d Sarek-data ]]
  then
    echo "$(tput setaf 1)Cloning Sarek-data repository$(tput sgr0)"
    git clone https://github.com/SciLifeLab/Sarek-data.git
  fi
  if [[ ! -d References ]]
  then
    echo "$(tput setaf 1)Building references$(tput sgr0)"
    nextflow run buildReferences.nf --refDir Sarek-data/reference --outDir References/$GENOME -profile $PROFILE --genome $GENOME --verbose
  fi
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
  run_wrapper --germline --sampleDir Sarek-data/testdata/tiny/normal
  clean_repo
fi

if [[ ALL,STEP =~ $TEST ]]
then
  run_wrapper --germline --sampleDir Sarek-data/testdata/tiny/normal
  run_wrapper --germline --step realign --noReports
  run_wrapper --germline --step recalibrate --noReports
  clean_repo
fi

if [[ ALL,GERMLINE =~ $TEST ]]
then
  run_wrapper --germline --sampleDir Sarek-data/testdata/tiny/normal --variantCalling --tools HaplotypeCaller
  clean_repo
fi

if [[ ALL,TOOLS =~ $TEST ]]
then
  run_wrapper --somatic --sample $SAMPLE --variantCalling  --tools FreeBayes,HaplotypeCaller,MuTect1,MuTect2
fi

if [[ ALL,MANTA =~ $TEST ]]
then
  run_wrapper --somatic --sample Sarek-data/testdata/tsv/tiny-manta.tsv --variantCalling --tools Manta --noReports
  run_wrapper --somatic --sample Sarek-data/testdata/tsv/tiny-manta.tsv --variantCalling --tools Manta,Strelka --noReports --strelkaBP
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
  run_wrapper --annotate --tools ${ANNOTATOR} --annotateVCF Sarek-data/testdata/vcf/Strelka_1234N_variants.vcf.gz --noReports
  run_wrapper --annotate --tools ${ANNOTATOR} --annotateVCF Sarek-data/testdata/vcf/Strelka_1234N_variants.vcf.gz,Sarek-data/testdata/vcf/Strelka_9876T_variants.vcf.gz
  clean_repo
fi

if [[ ALL,BUILDCONTAINERS =~ $TEST ]] && [[ $PROFILE == docker ]]
then
  ./scripts/do_all.sh --genome $GENOME
fi
