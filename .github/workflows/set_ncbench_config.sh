    #!/bin/bash

    for READS in 75 200; do
        yq --inplace '
            with(.variant-calls.nf-core-sarek-'"${PIPELINE_VERSION_NO_DOTS}"'-strelka-agilent-'"${READS}"'M.labels;
            .site = "nf-core" |
            .pipeline = "nf-core/sarek v'"${PIPELINE_VERSION}"'" |
            .trimming = "FastP v123" |
            .read-mapping = "bwa mem v123" |
            .base-quality-recalibration = "gatk4 v123" |
            .realignment = "none" |
            .variant-detection  = "strelka2 v123" |
            .genotyping = "none" |
            .reads = "'"${READS}"'M" ) |
            with(.variant-calls.nf-core-sarek-'"${PIPELINE_VERSION_NO_DOTS}"'-strelka-agilent-'"${READS}"'M.subcategory;
            . = "NA12878-agilent" ) |
            with(.variant-calls.nf-core-sarek-'"${PIPELINE_VERSION_NO_DOTS}"'-strelka-agilent-'"${READS}"'M.zenodo;
            .deposition = '"${DEPOSITION_ID}"'  |
            .filename = "get here proper file names" ) |
            with(.variant-calls.nf-core-sarek-'"${PIPELINE_VERSION_NO_DOTS}"'-strelka-agilent-'"${READS}"'M.benchmark;
            . = "giab-NA12878-agilent-'"${READS}"'M" ) |
            with(.variant-calls.nf-core-sarek-'"${PIPELINE_VERSION_NO_DOTS}"'-strelka-agilent-'"${READS}"'M.rename-contigs;
            . = "resources/rename-contigs/ucsc-to-ensembl.txt" )
            '  ncbench-workflow/config/config.yaml
    done



