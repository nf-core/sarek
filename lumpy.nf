#!/usr/bin/env nextflow
//This script runs tumor normal analysis using lumpy
//Author Jesper Eisfeldt
//Thanks to szilvajuhos and Pall for the nice "inspirational" code

params.sample = "lumpy"
// First checking the existence of the BAM file and its index for the tumor
params.tumor_bam = "tumor.bam" // override with --tumor_bam <SAMPLE>
params.tumor_bai = params.tumor_bam.replaceFirst(/.bam/,".bam.bai")
tumor_bam = file(params.tumor_bam)
if(!tumor_bam.exists()) exit 1, "Missing tumor file ${tumor_bam}; please specify --tumor_bam <TUMOR_BAM> --normal_bam <NORMAL_BAM> --sample <SAMPLE_ID>" 


tumor_bai = file(params.tumor_bai)
if(!tumor_bai.exists()) exit 1, "Missing tumor file index ${tumor_bai}" 

// Ditto for the normal
params.normal_bam = "normal.bam" // override with --normal_bam <SAMPLE>
normal_bam = file(params.normal_bam)
if(!normal_bam.exists()) exit 1, "Missing normal file ${normal_bam}; please specify --tumor_bam <TUMOR_BAM> --normal_bam <NORMAL_BAM> --sample <SAMPLE_ID>"
params.normal_bai = params.normal_bam.replaceFirst(/.bam/,".bam.bai")
normal_bai = file(params.normal_bai)
if(!normal_bai.exists()) exit 1, "Missing normal file index ${normal_bai}"

params.out = "$PWD"
//The uppmax module system does not include the split reads extraction script
params.splitreads_extract="/sw/apps/bioinfo/LUMPY/0.2.12/milou/scripts/extractSplitReads_BwaMem"

process lumpy{

    module 'bioinfo-tools'
    module 'samtools'
    module 'LUMPY/0.2.12'

    input:

    file tumor_bam
    file normal_bam

    file tumor_bai
    file normal_bai

    output:
       file "${params.sample}.tumor_normal.lumpySV.vcf" into lumpy_vcf

    """
    samtools view -b -F 1294 ${params.tumor_bam} > ${params.sample}_tumor.D.unsorted.bam
    samtools view -b -F 1294 ${params.normal_bam} > ${params.sample}_normal.D.unsorted.bam

    samtools sort ${params.sample}_tumor.D.unsorted.bam ${params.sample}_tumor.D
    samtools sort ${params.sample}_normal.D.unsorted.bam ${params.sample}_normal.D

    rm ${params.sample}_tumor.D.unsorted.bam ${params.sample}_normal.D.unsorted.bam

    samtools view -h ${params.tumor_bam} | ${params.splitreads_extract} -i stdin | samtools view -Sb - > ${params.sample}_tumor.S.unsorted.bam
    samtools view -h ${params.normal_bam} | ${params.splitreads_extract} -i stdin | samtools view -Sb - > ${params.sample}_normal.S.unsorted.bam

    samtools sort ${params.sample}_tumor.S.unsorted.bam ${params.sample}_tumor.S
    samtools sort ${params.sample}_normal.S.unsorted.bam ${params.sample}_normal.S
    rm ${params.sample}_tumor.S.unsorted.bam ${params.sample}_normal.S.unsorted.bam

    lumpyexpress -B ${params.normal_bam},${params.tumor_bam} -S ${params.sample}_normal.S.bam,${params.sample}_tumor.S.bam -D ${params.sample}_normal.D.bam,${params.sample}_tumor.D.bam -o ${params.sample}.tumor_normal.lumpySV.vcf
    rm ${params.sample}_normal.S.bam ${params.sample}_tumor.S.bam ${params.sample}_normal.D.bam ${params.sample}_tumor.D.bam
    """
}
