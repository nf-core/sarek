#!/bin/bash
echo "Starting Nextflow... Command:"
echo "./nextflow run buildReferences.nf -profile testing --storeDirectory smallGRCh37"
echo "-----"
./nextflow run buildReferences.nf -profile testing --storeDirectory smallGRCh37
echo "Cleaning up docker images:"
echo "docker rmi -f $(docker images -q)"
echo "-----"
docker rmi -f $(docker images -q)
echo "Starting Nextflow... Command:"
echo "./nextflow run . -profile testing --test --step preprocessing --tools MultiQC"
echo "-----"
./nextflow run . -profile testing --test --step preprocessing --tools MultiQC
echo "Cleaning up docker images:"
echo "docker rmi -f $(docker images -q)"
echo "-----"
docker rmi -f $(docker images -q)
echo "Starting Nextflow... Command:"
echo "./nextflow run . -profile testing --test --step realign"
echo "-----"
./nextflow run . -profile testing --test --step realign
echo "Cleaning up docker images:"
echo "docker rmi -f $(docker images -q)"
echo "-----"
docker rmi -f $(docker images -q)
echo "Starting Nextflow... Command:"
echo "./nextflow run . -profile testing --test --step recalibrate"
echo "-----"
./nextflow run . -profile testing --test --step recalibrate
echo "Cleaning up docker images:"
echo "docker rmi -f $(docker images -q)"
echo "-----"
docker rmi -f $(docker images -q)
echo "Starting Nextflow... Command:"
echo "./nextflow run . -profile testing --test --step skipPreprocessing --tools HaplotypeCaller,MultiQC,MuTect1,MuTect2,Strelka,snpEff,VEP"
echo "-----"
./nextflow run . -profile testing --test --step skipPreprocessing --tools HaplotypeCaller,MultiQC,MuTect1,MuTect2,Strelka,snpEff,VEP
