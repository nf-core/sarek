#!/bin/bash
echo "Starting Nextflow... Command:"
echo "nextflow run buildReferences.nf -profile docker --storeDirectory smallGRCh37"
echo "-----"
nextflow run buildReferences.nf -profile docker --storeDirectory smallGRCh37
echo "Starting Nextflow... Command:"
echo "nextflow run . -profile docker --test --step preprocessing --tools MultiQC"
echo "-----"
nextflow run . -profile docker --test --step preprocessing --tools MultiQC
echo "Starting Nextflow... Command:"
echo "nextflow run . -profile docker --test --step realign --tools MultiQC "
echo "-----"
nextflow run . -profile docker --test --step realign --tools MultiQC
echo "Starting Nextflow... Command:"
echo "nextflow run . -profile docker --test --step recalibrate --tools MultiQC "
echo "-----"
nextflow run . -profile docker --test --step recalibrate --tools MultiQC
echo "Starting Nextflow... Command:"
echo "nextflow run . -profile docker --test --step skipPreprocessing --tools HaplotypeCaller,MuTect1,MuTect2 "
echo "-----"
nextflow run . -profile docker --test --step skipPreprocessing --tools HaplotypeCaller,MuTect1,MuTect2
