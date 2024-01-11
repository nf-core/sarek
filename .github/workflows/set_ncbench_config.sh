#!/bin/bash
#        .trimming = "FastP v${FASTP_VERSION}" |
#        .read-mapping = "bwa mem v${BWA_VERSION}" |
#        .base-quality-recalibration = "gatk4 v${GATK_VERSION}" |
#        .realignment = "none" |
#        .variant-detection  = "strelka2 v${STRELKA_VERSION}" |
#        .genotyping = "none" |
for READS in 75 200; do
    yq --inplace "
        with(.variant-calls.nf-core-sarek-$PIPELINE_VERSION-strelka-agilent-$READSM.labels;
        .site = "nf-core" |
        .pipeline = "nf-core/sarek v$PIPELINE_VERSION" |
        .reads = "$READSM" ) |
        with(.variant-calls.nf-core-sarek-$PIPELINE_VERSION-strelka-agilent-$READSM.subcategory;
        . = "NA12878-agilent" ) |
        with(.variant-calls.nf-core-sarek-$PIPELINE_VERSION-strelka-agilent-$READSM.zenodo;
        .deposition = $DEPOSITION_ID  |
        .filename = "TODO get here proper file names" ) |
        with(.variant-calls.nf-core-sarek-$PIPELINE_VERSION-strelka-agilent-$READSM.benchmark;
        . = "giab-NA12878-agilent-$READSM" ) |
        with(.variant-calls.nf-core-sarek-$PIPELINE_VERSION-strelka-agilent-$READSM.rename-contigs;
        . = "resources/rename-contigs/ucsc-to-ensembl.txt" )
        " ncbench-workflow/config/config.yaml
done



