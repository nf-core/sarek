#!/bin/bash
echo "Starting Nextflow... Command:"
echo "./nextflow run buildReferences.nf -profile testing"
echo "-----"
./nextflow run buildReferences.nf -profile testing
echo "Cleaning up docker images:"
echo "docker rmi -f maxulysse/igvtools:1.1"
echo "-----"
docker rmi -f maxulysse/igvtools:1.1
echo "Starting Nextflow... Command:"
echo "./nextflow run . -profile testing --test --step preprocessing --tools MultiQC"
echo "-----"
./nextflow run . -profile testing --test --step preprocessing --tools MultiQC
echo "Cleaning up docker images:"
echo "docker rmi -f maxulysse/mapreads:1.1 maxulysse/samtools:1.1 maxulysse/picard:1.1 maxulysse/fastqc:1.1"
echo "-----"
docker rmi -f maxulysse/mapreads:1.1 maxulysse/samtools:1.1 maxulysse/picard:1.1 maxulysse/fastqc:1.1
echo "Starting Nextflow... Command:"
echo "./nextflow run . -profile testing --test --step realign"
echo "-----"
./nextflow run . -profile testing --test --step realign
echo "Starting Nextflow... Command:"
echo "./nextflow run . -profile testing --test --step recalibrate"
echo "-----"
./nextflow run . -profile testing --test --step recalibrate
echo "Starting Nextflow... Command:"
echo "./nextflow run . -profile testing --test --step skipPreprocessing --tools HaplotypeCaller,MultiQC,MuTect1,MuTect2,Strelka,snpEff,VEP"
echo "-----"
./nextflow run . -profile testing --test --step skipPreprocessing --tools HaplotypeCaller,MultiQC,MuTect1,MuTect2,Strelka,snpEff,VEP
