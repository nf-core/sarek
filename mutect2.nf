#!/usr/bin/env nextflow

params.tumor_bam =      "/proj/a2014205/nobackup/szilva/TCGA/S1/downsampling/G15511/recalibrated/G15511.normal__0.recal.bam"
params.tumor_bam_idx =  "/proj/a2014205/nobackup/szilva/TCGA/S1/downsampling/G15511/recalibrated/G15511.normal__0.recal.bai"
params.normal_bam =     "/proj/a2014205/nobackup/szilva/TCGA/S1/downsampling/G15511/recalibrated/G15511.tumor__1.recal.bam"
params.normal_bam_idx = "/proj/a2014205/nobackup/szilva/TCGA/S1/downsampling/G15511/recalibrated/G15511.tumor__1.recal.bai"

tumor_bam = file(params.tumor_bam)
tumor_bam_idx = file(params.tumor_bam_idx)
normal_bam = file(params.normal_bam)
normal_bam_idx = file(params.normal_bam_idx)

// define intervals file by --intervals
intervalsFile = file(params.intervals)
intervals = Channel
    .from(intervalsFile.readLines())

process runIntervals {

    module 'bioinfo-tools'
    module 'java/sun_jdk1.8.0_92'

    threads 8
    memory { 32.GB * task.attempt }
    time { 16.h * task.attempt }
    errorStrategy { task.exitStatus == 143 ? 'retry' : 'terminate' }
    maxRetries 3
    maxErrors '-1'

    input:
    file tumor_bam
    file tumor_bam_idx
    file normal_bam 
    file normal_bam_idx
    set genomicInterval from intervals

    output:
    file '*.vcf' into finalsVCFs

    // we are using MuTect2 shipped in GATK v3.6

    """
    java -Xmx${task.memory.toGiga()}g -jar /home/szilva/dev/GATK/GenomeAnalysisTK.jar \
    -T MuTect2 \
    -nct ${task.threads} \
    -R /sw/data/uppnex/ToolBox/ReferenceAssemblies/hg38make/bundle/2.8/b37/human_g1k_v37_decoy.fasta \
    -I:tumor ${tumor_bam}\
    -I:normal ${normal_bam} \
    --dbsnp /sw/data/uppnex/ToolBox/ReferenceAssemblies/hg38make/bundle/2.8/b37/dbsnp_138.b37.vcf \
    --cosmic /sw/data/uppnex/ToolBox/ReferenceAssemblies/hg38make/bundle/2.8/b37/b37_cosmic_v54_120711.vcf \
    -L \"$genomicInterval\" \
    -o ${genomicInterval}.vcf
    """
}
