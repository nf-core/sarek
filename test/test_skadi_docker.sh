#!/bin/bash
echo "Starting Nextflow... Command:"
echo "nextflow run . -profile docker_skadi --test --steps preprocessing --tools MultiQC"
echo "-----"
nextflow run . -profile docker_skadi --test --steps preprocessing --tools MultiQC
echo "Starting Nextflow... Command:"
echo "nextflow run . -profile docker_skadi --test --steps realign --tools MultiQC"
echo "-----"
nextflow run . -profile docker_skadi --test --steps realign --tools MultiQC
echo "Starting Nextflow... Command:"
echo "nextflow run . -profile docker_skadi --test --steps recalibrate --tools MultiQC"
echo "-----"
nextflow run . -profile docker_skadi --test --steps recalibrate --tools MultiQC
echo "Starting Nextflow... Command:"
echo "nextflow run . -profile docker_skadi --test --steps skipPreprocessing --tools HaplotypeCaller,MuTect1,MuTect2"
echo "-----"
nextflow run . -profile docker_skadi --test --steps skipPreprocessing --tools HaplotypeCaller,MuTect1,MuTect2
