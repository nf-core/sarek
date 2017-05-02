#!/bin/bash
echo "Starting Nextflow... Command:"
echo "./nextflow run buildReferences.nf -profile testing --download"
echo "-----"
./nextflow run buildReferences.nf -profile testing --download
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
echo "./nextflow run . -profile testing --step realign"
echo "-----"
./nextflow run . -profile testing --step realign
echo "Starting Nextflow... Command:"
echo "./nextflow run . -profile testing --step recalibrate"
echo "-----"
./nextflow run . -profile testing --step recalibrate
echo "Starting Nextflow... Command:"
echo "./nextflow run . -profile testing --step skipPreprocessing --tools FreeBayes,HaplotypeCaller,MultiQC,MuTect1,MuTect2,Strelka,VarDict"
echo "-----"
./nextflow run . -profile testing --step skipPreprocessing --tools FreeBayes,HaplotypeCaller,MultiQC,MuTect1,MuTect2,Strelka,VarDict
echo "Cleaning up docker images:"
echo "docker rmi -f $(docker images -q)"
echo "-----"
docker rmi -f $(docker images -q)
echo "Starting Nextflow... Command:"
echo "./nextflow run . -profile testing --step skipPreprocessing --tools MuTect2,snpEff"
echo "-----"
./nextflow run . -profile testing --step skipPreprocessing --tools MuTect2,snpEff
echo "Starting Nextflow... Command:"
echo "./nextflow run . -profile testing --step annotate --tools snpEff --annotateTools MuTect2"
echo "-----"
./nextflow run . -profile testing --step annotate --tools snpEff --annotateTools MuTect2
# echo "Starting Nextflow... Command:"
# echo "./nextflow run . -profile testing --step annotate --tools snpEff --annotateVCF VariantCalling/MuTect2/mutect2_9876T_vs_1234N.vcf.gz,VariantCalling/MuTect2/mutect2_9877R_vs_1234N.vcf.gz"
# echo "-----"
# ./nextflow run . -profile testing --step annotate --tools snpEff --annotateVCF VariantCalling/MuTect2/mutect2_9876T_vs_1234N.vcf.gz,VariantCalling/MuTect2/mutect2_9877R_vs_1234N.vcf.gz
