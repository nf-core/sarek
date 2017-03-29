#!/bin/bash
echo "Starting Nextflow... Command:"
echo "./nextflow run buildReferences.nf -profile testing --storeDirectory smallGRCh37"
echo "-----"
./nextflow run buildReferences.nf -profile testing --storeDirectory smallGRCh37
echo "Starting Nextflow... Command:"
echo "./nextflow run . -profile testing --test --step preprocessing --tools MultiQC"
echo "-----"
./nextflow run . -profile testing --test --step preprocessing --tools MultiQC
echo "Starting Nextflow... Command:"
echo "./nextflow run . -profile testing --test --step realign"
echo "-----"
./nextflow run . -profile testing --test --step realign
echo "Starting Nextflow... Command:"
echo "./nextflow run . -profile testing --test --step recalibrate"
echo "-----"
./nextflow run . -profile testing --test --step recalibrate
echo "Starting Nextflow... Command:"
echo "./nextflow run . -profile testing --test --step skipPreprocessing --tools HaplotypeCaller,MultiQC,MuTect1,MuTect2,Strelka,snpEff"
echo "-----"
./nextflow run . -profile testing --test --step skipPreprocessing --tools HaplotypeCaller,MultiQC,MuTect1,MuTect2,Strelka,snpEff
