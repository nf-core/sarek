#!/bin/bash

declare -A variant_callers=(
    ["deepvariant"]="deepvariant/NA12878_%sM/NA12878_%sM.deepvariant.vcf.gz"
    ["freebayes"]="freebayes/NA12878_%sM/NA12878_%sM.freebayes.vcf.gz"
    ["haplotypecaller"]="haplotypecaller/NA12878_%sM/NA12878_%sM.haplotypecaller.filtered.vcf.gz"
    ["strelka2"]="strelka/NA12878_%sM/NA12878_%sM.strelka.variants.vcf.gz"
)


declare -A variant_versions=(
    ["deepvariant"]="${DEEPVARIANT_VERSION}"
    ["freebayes"]="${FREEBAYES_VERSION}"
    ["haplotypecaller"]="${HAPLOTYPECALLER_VERSION}"
    ["strelka2"]="${STRELKA_VERSION}"
)

for READS in 75 200; do
    for variant_caller in "${!variant_callers[@]}"; do
        filename=$(printf "${variant_callers[$variant_caller]}" $READS $READS)
        yq --inplace '
            with(.variant-calls.nf-core-sarek-'"${PIPELINE_VERSION_NO_DOTS}"'-'"${variant_caller}"'-agilent-'"${READS}"'M.labels;
            .site = "nf-core" |
            .pipeline = "nf-core/sarek v'"${PIPELINE_VERSION}"'" |
            .trimming = "FastP v'"${FASTP_VERSION}"'" |
            .read-mapping = "bwa mem v'"${BWA_VERSION}"'" |
            .base-quality-recalibration = "gatk4 v'"${BQSR_VERSION}"'" |
            .realignment = "none" |
            .variant-detection  = "'${variant_caller}' v'"${variant_versions[$variant_caller]}"'" |
            .genotyping = "none" |
            .reads = "'"${READS}"'M" ) |
            with(.variant-calls.nf-core-sarek-'"${PIPELINE_VERSION_NO_DOTS}"'-'"${variant_caller}"'-agilent-'"${READS}"'M.subcategory;
            . = "NA12878-agilent" ) |
            with(.variant-calls.nf-core-sarek-'"${PIPELINE_VERSION_NO_DOTS}"'-'"${variant_caller}"'-agilent-'"${READS}"'M.zenodo;
            .deposition = '"${DEPOSITION_ID}"'  |
            .filename= "'"${filename}"'" ) |
            with(.variant-calls.nf-core-sarek-'"${PIPELINE_VERSION_NO_DOTS}"'-'"${variant_caller}"'-agilent-'"${READS}"'M.benchmark;
            . = "giab-NA12878-agilent-'"${READS}"'M" ) |
            with(.variant-calls.nf-core-sarek-'"${PIPELINE_VERSION_NO_DOTS}"'-'"${variant_caller}"'-agilent-'"${READS}"'M.rename-contigs;
            . = "resources/rename-contigs/ucsc-to-ensembl.txt" )
            '  ncbench-workflow/config/config.yaml
    done
done

