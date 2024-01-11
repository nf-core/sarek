#!/bin/bash

        #with(.variant-calls.nf-core-sarek-$PIPELINE_VERSION-strelka-agilent-${READS}M.zenodo;
        #.deposition = $DEPOSITION_ID  |
        #.filename = "TODO get here proper file names" ) |

        # .trimming = "FastP v$FASTP_VERSION" |
        # .read-mapping = "bwa mem v$BWA_VERSION" |
        # .base-quality-recalibration = "gatk4 v$GATK_VERSION" |
        # .realignment = "none" |
        # .variant-detection  = "strelka2 v$STRELKA_VERSION" |

        #with(.variant-calls.nf-core-sarek-$PIPELINE_VERSION-strelka-agilent-${READS}M.subcategory;
        #. = "NA12878-agilent" ) |
        #with(.variant-calls.nf-core-sarek-$PIPELINE_VERSION-strelka-agilent-${READS}M.benchmark;
        #. = "giab-NA12878-agilent-${READS}M" ) |

        #with(.variant-calls.nf-core-sarek-$PIPELINE_VERSION-strelka-agilent-${READS}M.rename-contigs;
        #. = "resources/rename-contigs/ucsc-to-ensembl.txt" )

ls ncbench-workflow

for READS in 75 200; do
    yq --inplace "
        with(.variant-calls.nf-core-sarek-$PIPELINE_VERSION-strelka-agilent-${READS}M.labels;
        .site = "nf-core" |
        .pipeline = "nf-core/sarek v$PIPELINE_VERSION" |
        .genotyping = "none" |
        .reads = "${READS}M" )
        " ncbench-workflow/config/config.yaml
done



