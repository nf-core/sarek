#!/bin/bash
set -xeuo pipefail

GENOME="smallGRCh37"
PROFILE="singularityTest"
TEST="ALL"
TRAVIS=${TRAVIS:-false}
SAMPLE="data/tsv/tiny.tsv"

while [[ $# -gt 0 ]]
do
  key="$1"
  case $key in
    -g|--genome)
    GENOME="$2"
    shift
    ;;
    -p|--profile)
    PROFILE="$2"
    shift
    ;;
    -s|--sample)
    SAMPLE="$2"
    shift
    ;;
    -t|--test)
    TEST="$2"
    shift
    ;;
    *) # unknown option
    ;;
  esac
  shift
done

function nf_test() {
  echo "$(tput setaf 1)nextflow run $@ -profile $PROFILE --genome $GENOME -resume --verbose$(tput sgr0)"
  nextflow run "$@" -profile "$PROFILE" --genome $GENOME -resume --verbose
}

# Build references only for smallGRCh37
if [[ "$GENOME" == "smallGRCh37" ]] && [[ "$TEST" != "BUILDCONTAINERS" ]]
then
  nf_test buildReferences.nf --download

  # Remove images only on TRAVIS
  if [[ "$PROFILE" == "dockerTest" ]] && [[ "$TRAVIS" == true ]]
  then
    docker rmi -f maxulysse/igvtools:1.1
  elif [[ "$PROFILE" == singularityTest ]] && [[ "$TRAVIS" == true ]]
  then
    rm -rf work/singularity/igvtools-1.1.img
  fi
fi

if [[ "$TEST" = "MAPPING" ]] || [[ "$TEST" = "ALL" ]]
then
  nf_test . --step preprocessing --sample $SAMPLE
fi

if [[ "$TEST" = "REALIGN" ]] || [[ "$TEST" = "ALL" ]]
then
  nf_test . --step preprocessing --sample $SAMPLE
  nf_test . --step realign --noReports
  nf_test . --step realign --tools HaplotypeCaller
  nf_test . --step realign --tools HaplotypeCaller --noReports --noGVCF
fi

if [[ "$TEST" = "RECALIBRATE" ]] || [[ "$TEST" = "ALL" ]]
then
  nf_test . --step preprocessing --sample $SAMPLE
  nf_test . --step recalibrate --noReports
  nf_test . --step recalibrate --tools FreeBayes,HaplotypeCaller,MuTect1,MuTect2,Strelka
  # Test whether restarting from an already recalibrated BAM works
  nf_test . --step skipPreprocessing --tools Strelka --noReports
fi

if [[ "$TEST" = "ANNOTATEVEP" ]] || [[ "$TEST" = "ALL" ]]
then
  nf_test . --step preprocessing --sample data/tsv/tiny-single-manta.tsv --tools Manta
  nf_test . --step preprocessing --sample data/tsv/tiny-manta.tsv --tools Manta
  nf_test . --step preprocessing --sample $SAMPLE --tools MuTect2,Strelka

  # Remove images only on TRAVIS
  if [[ "$PROFILE" == "dockerTest" ]] && [[ "$TRAVIS" == true ]]
  then
    docker rmi -f maxulysse/fastqc:1.1 maxulysse/gatk:1.0 maxulysse/gatk:1.1 maxulysse/mapreads:1.1 maxulysse/picard:1.1 maxulysse/runmanta:dev maxulysse/samtools:1.1 maxulysse/strelka:dev
  elif [[ "$PROFILE" == "singularityTest" ]] && [[ "$TRAVIS" == true ]]
  then
    rm -rf work/singularity/fastqc-1.1.img work/singularity/gatk-1.0.img work/singularity/gatk-1.1.img work/singularity/mapreads-1.1.img work/singularity/picard-1.1.img work/singularity/runmanta-dev.img work/singularity/samtools-1.1.img work/singularity/strelka-dev.img
  fi
  nf_test . --step annotate --tools VEP --annotateTools Manta,Strelka
  nf_test . --step annotate --tools VEP --annotateVCF VariantCalling/Manta/Manta_9876T_vs_1234N.diploidSV.vcf.gz,VariantCalling/Manta/Manta_9876T_vs_1234N.somaticSV.vcf.gz --noReports
  nf_test . --step annotate --tools VEP --annotateVCF VariantCalling/Manta/Manta_9876T_vs_1234N.diploidSV.vcf.gz --noReports
fi

if [[ "$TEST" = "ANNOTATESNPEFF" ]] || [[ "$TEST" = "ALL" ]]
then
  nf_test . --step preprocessing --sample data/tsv/tiny-single-manta.tsv --tools Manta
  nf_test . --step preprocessing --sample data/tsv/tiny-manta.tsv --tools Manta
  nf_test . --step preprocessing --sample $SAMPLE --tools MuTect2,Strelka

  # Remove images only on TRAVIS
  if [[ "$PROFILE" == "dockerTest" ]] && [[ "$TRAVIS" == true ]]
  then
    docker rmi -f maxulysse/fastqc:1.1 maxulysse/gatk:1.0 maxulysse/gatk:1.1 maxulysse/mapreads:1.1 maxulysse/picard:1.1 maxulysse/runmanta:dev maxulysse/samtools:1.1 maxulysse/strelka:dev
  elif [[ "$PROFILE" == "singularityTest" ]] && [[ "$TRAVIS" == true ]]
  then
    rm -rf work/singularity/fastqc-1.1.img work/singularity/gatk-1.0.img work/singularity/gatk-1.1.img work/singularity/mapreads-1.1.img work/singularity/picard-1.1.img work/singularity/runmanta-dev.img work/singularity/samtools-1.1.img work/singularity/strelka-dev.img
  fi
  nf_test . --step annotate --tools snpEff --annotateTools Manta,Strelka
  nf_test . --step annotate --tools snpEff --annotateVCF VariantCalling/Manta/Manta_9876T_vs_1234N.diploidSV.vcf.gz,VariantCalling/Manta/Manta_9876T_vs_1234N.somaticSV.vcf.gz --noReports
  nf_test . --step annotate --tools snpEff --annotateVCF VariantCalling/Manta/Manta_9876T_vs_1234N.diploidSV.vcf.gz --noReports
fi

if [[ "$TEST" = "BUILDCONTAINERS" ]] || [[ "$TEST" = "ALL" ]]
then
  nf_test buildContainers.nf --docker --containers bcftools,concatvcf,fastqc,gatk,htslib,igvtools,mapreads,multiqc,mutect1,picard,qualimap,runallelecount,runascat,runconvertallelecounts,runmanta,samtools,strelka,snpeff,vep
fi
