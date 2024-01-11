
#!/bin/bash

for reads in 75 200; do
    yq --inplace '
        with(.variant-calls.nf-core-sarek-${{ env.PIPELINE_VERSION }}-strelka-agilent-${reads}M.labels;
        .site = "nf-core" |
        .pipeline = "nf-core/sarek v${{ env.PIPELINE_VERSION }}" |
        .trimming = "FastP v${{ env.FASTP_VERSION }}" |
        .read-mapping = "bwa mem v${{ env.BWA_VERSION }}" |
        .base-quality-recalibration = "gatk4 v${{ env.GATK_VERSION }}" |
        .realignment = "none" |
        .variant-detection  = "strelka2 v${{ env.STRELKA_VERSION }}" |
        .genotyping = "none" |
        .reads = "${reads}M" ) |
        with(.variant-calls.nf-core-sarek-${{ env.PIPELINE_VERSION }}-strelka-agilent-${reads}M.subcategory;
        . = "NA12878-agilent" ) |
        with(.variant-calls.nf-core-sarek-${{ env.PIPELINE_VERSION }}-strelka-agilent-${reads}M.zenodo;
        .deposition = $DEPOSITION_ID  |
        .filename = "TODO get here proper file names" ) |
        with(.variant-calls.nf-core-sarek-${{ env.PIPELINE_VERSION }}-strelka-agilent-${reads}M.benchmark;
        . = "giab-NA12878-agilent-${reads}M" ) |
        with(.variant-calls.nf-core-sarek-${{ env.PIPELINE_VERSION }}-strelka-agilent-${reads}M.rename-contigs;
        . = "resources/rename-contigs/ucsc-to-ensembl.txt" )
        ' ncbench-workflow/config/config.yaml
done



