#!/bin/bash
set -xeuo pipefail

BUILD=false
KEEP=false
GENOME=smallGRCh37
PROFILE=singularity
SAMPLE=Sarek-data/testdata/tsv/tiny.tsv
TEST=ALL
TRAVIS=${TRAVIS:-false}
CPUS=2

TMPDIR=`pwd`/tmp
mkdir -p $TMPDIR
export NXF_SINGULARITY_CACHEDIR=$TMPDIR
export NXF_TEMP=$TMPDIR
export SINGULARITY_CACHEDIR=$TMPDIR
export SINGULARITY_TMPDIR=$TMPDIR

# remove Reference directory
rm -rf References

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
    -c|--cpus)
    CPUS=$2
    shift # past value
    ;;
    *) # unknown option
    shift # past argument
    ;;
  esac
done

function run_wrapper() {
  ./scripts/wrapper.sh $@ --profile $PROFILE --genome $GENOME --genomeBase $PWD/References/$GENOME --verbose --cpus ${CPUS}
}

function clean_repo() {
  if [[ $TRAVIS == false ]] && [[ $KEEP == false ]]
  then
    echo "$(tput setaf 1)Cleaning directory$(tput sgr0)"
    rm -rf work .nextflow* Annotation Preprocessing Reports Results VariantCalling
  fi
}

# Build references only for smallGRCh37
if [[ $GENOME == smallGRCh37 ]] && [[ $TEST != BUILDCONTAINERS ]] && [[ BUILD ]]
then
  if [[ -z "$(ls -A Sarek-data)" ]]
  then
    git submodule init
    git submodule update
  fi
  if [[ ! -d References ]]
  then
    echo "$(tput setaf 1)Building references$(tput sgr0)"
    nextflow run buildReferences.nf --refDir Sarek-data/reference --outDir References/$GENOME -profile $PROFILE --genome $GENOME --verbose
  fi
fi


if [[ ALL,GERMLINE =~ $TEST ]]
then
	# Added Strelka to germline test (no Strelka best practices test for this small data) and not asking for reports
	run_wrapper --germline --sampleDir Sarek-data/testdata/tiny/normal --variantCalling --tools HaplotypeCaller,Strelka --noReports
	run_wrapper --germline --sampleDir Sarek-data/testdata/tiny/normal --variantCalling --tools HaplotypeCaller,Strelka --bed Sarek-data/testdata/target.bed --noReports
	run_wrapper --germline --step recalibrate --noReports
	clean_repo
fi

if [[ ALL,SOMATIC =~ $TEST ]]
then
	run_wrapper --somatic --sample Sarek-data/testdata/tsv/tiny-manta.tsv --variantCalling --tools FreeBayes,HaplotypeCaller,Manta,Mutect2 --noReports
	run_wrapper --somatic --sample Sarek-data/testdata/tsv/tiny-manta.tsv --variantCalling --tools Manta,Strelka --noReports --strelkaBP
	clean_repo
fi

if [[ ALL,TARGETED =~ $TEST ]]
then
	run_wrapper --somatic --sample Sarek-data/testdata/tsv/tiny-manta.tsv --variantCalling --tools Manta,Strelka --noReports --bed Sarek-data/testdata/target.bed
	clean_repo
fi

if [[ ALL,ANNOTATEALL,ANNOTATESNPEFF,ANNOTATEVEP =~ $TEST ]]
then
  if [[ $TEST = ANNOTATESNPEFF ]]
  then
    ANNOTATOR=snpEFF
  elif [[ $TEST = ANNOTATEVEP ]]
  then
    ANNOTATOR=VEP
  elif [[ ALL,ANNOTATEALL =~ $TEST ]]
  then
    ANNOTATOR=merge,snpEFF,VEP
  fi
  run_wrapper --annotate --tools ${ANNOTATOR} --annotateVCF Sarek-data/testdata/vcf/Strelka_1234N_variants.vcf.gz
  clean_repo
fi

if [[ MULTIPLE =~ $TEST ]]
then
  run_wrapper --somatic --sample Sarek-data/testdata/tsv/tiny-multiple.tsv --variantCalling --tools FreeBayes,HaplotypeCaller,Manta,Mutect2 --noReports
	run_wrapper --somatic --sample Sarek-data/testdata/tsv/tiny-multiple.tsv --variantCalling --tools Manta,Strelka --noReports --strelkaBP
fi

if [[ BUILDCONTAINERS =~ $TEST ]] && [[ $PROFILE == docker ]]
then
  ./scripts/do_all.sh --genome $GENOME
fi
