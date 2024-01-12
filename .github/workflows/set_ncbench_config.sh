#!/bin/bash

for READS in 75 200; do
    yq --inplace '
            with(.variant-calls.nf-core-sarek-${PIPELINE_VERSION_NO_DOTS}-strelka-agilent-75M.labels;
            .site = "nf-core" |
            .pipeline = "nf-core/sarek v${PIPELINE_VERSION}" |
            .trimming = "FastP v${FASTP_VERSION }" |
            .read-mapping = "bwa mem v${ BWA_VERSION }" |
            .base-quality-recalibration = "gatk4 v${ GATK_VERSION }" |ß
            .realignment = "none" |
            .variant-detection  = "strelka2 v${ STRELKA_VERSION }" |
            .genotyping = "none" |
            .reads = "75M" ) |
            with(.variant-calls.nf-core-sarek-${ PIPELINE_VERSION_NO_DOTS }-strelka-agilent-75M.subcategory;
            . = "NA12878-agilent" ) |
            with(.variant-calls.nf-core-sarek-${ PIPELINE_VERSION_NO_DOTS }-strelka-agilent-75M.zenodo;
            .deposition = ${ DEPOSITION_ID }  |
            .filename = "get here proper file names" ) |
            with(.variant-calls.nf-core-sarek-${ PIPELINE_VERSION_NO_DOTS }-strelka-agilent-75M.benchmark;
            . = "giab-NA12878-agilent-75M" ) |
            with(.variant-calls.nf-core-sarek-${ PIPELINE_VERSION_NO_DOTS }-strelka-agilent-75M.rename-contigs;
            . = "resources/rename-contigs/ucsc-to-ensembl.txt" )
            '  ncbench-workflow/config/config.yaml
done



