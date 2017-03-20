#!/bin/bash
echo "Starting Nextflow... Command:"
echo "nextflow run . -profile docker_skadi --test --step preprocessing --tools MultiQC"
echo "-----"
nextflow run . -profile docker_skadi --test --step preprocessing --tools MultiQC
echo "Starting Nextflow... Command:"
echo "nextflow run . -profile docker_skadi --test --step realign --tools MultiQC"
echo "-----"
nextflow run . -profile docker_skadi --test --step realign --tools MultiQC
echo "Starting Nextflow... Command:"
echo "nextflow run . -profile docker_skadi --test --step recalibrate --tools MultiQC"
echo "-----"
nextflow run . -profile docker_skadi --test --step recalibrate --tools MultiQC
echo "Starting Nextflow... Command:"
echo "nextflow run . -profile docker_skadi --test --step skipPreprocessing --tools HaplotypeCaller,MuTect1,MuTect2"
echo "-----"
nextflow run . -profile docker_skadi --test --step skipPreprocessing --tools HaplotypeCaller,MuTect1,MuTect2
