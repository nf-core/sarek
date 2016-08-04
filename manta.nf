#!/usr/bin/env nextflow
//This script runs tumor normal analysis using manta
//Author Jesper Eisfeldt
//Thanks to szilvajuhos and Pall for the nice "inspirational" code

params.sample = "M.alfredi"
// First checking the existence of the BAM file and its index for the tumor
params.tumor_bam = "tumor.bam" // override with --tumor_bam <SAMPLE>
params.tumor_bai = params.tumor_bam.replaceFirst(/.bam/,".bam.bai")
tumor_bam = file(params.tumor_bam)
if(!tumor_bam.exists()) exit 1, "Missing tumor file ${tumor_bam}; please specify --tumor_bam <TUMOR_BAM> --normal_bam <NORMAL_BAM> --genome <REFERENCE_FASTA> --sample <SAMPLE_ID>" 


tumor_bai = file(params.tumor_bai)
if(!tumor_bai.exists()) exit 1, "Missing tumor file index ${tumor_bai}" 

// Ditto for the normal
params.normal_bam = "normal.bam" // override with --normal_bam <SAMPLE>
normal_bam = file(params.normal_bam)
if(!normal_bam.exists()) exit 1, "Missing normal file ${normal_bam}; please specify --tumor_bam <TUMOR_BAM> --normal_bam <NORMAL_BAM> --genome <REFERENCE_FASTA> --sample <SAMPLE_ID>"
params.normal_bai = params.normal_bam.replaceFirst(/.bam/,".bam.bai")
normal_bai = file(params.normal_bai)
if(!normal_bai.exists()) exit 1, "Missing normal file index ${normal_bai}"


params.genome = "/sw/data/uppnex/ToolBox/ReferenceAssemblies/hg38make/bundle/2.8/b37/human_g1k_v37_decoy.fasta"
genome_file= file(params.genome)
if(!genome_file.exists()) exit 1, "Missing reference file ${reference_fa}; please specify --tumor_bam <TUMOR_BAM> --normal_bam <NORMAL_BAME> --genome <REFERENCE_FASTA> --sample <SAMPLE_ID>"

params.genomeidx = "${params.genome}.fai"
genome_index = file(params.genomeidx)
if(!genome_index.exists()) exit 1, "Error: the fasta file must be indexed"



process manta{

    module 'bioinfo-tools'
    module 'manta'

    cpus 8

    input:
    file genome_file
    file genome_index

    file tumor_bam
    file normal_bam

    file tumor_bai
    file normal_bai

    output:
       file "${params.sample}.somaticSV.vcf" into manta_somatic_vcf
       file "${params.sample}.diploidSV.vcf" into manta_diploid_vcf
       file "${params.sample}.candidateSV.vcf" into manta_candidate_vcf
       file "${params.sample}.candidateSmallIndels.vcf" into manta_indel_vcf

    """
    configManta.py --normalBam ${params.normal_bam} --tumorBam ${params.tumor_bam} --reference ${params.genome} --runDir ${params.sample}_manta_dir
    python ${params.sample}_manta_dir/runWorkflow.py -m local -j 8
    gunzip -c ${params.sample}_manta_dir/results/variants/somaticSV.vcf.gz > ${params.sample}.somaticSV.vcf
    gunzip -c ${params.sample}_manta_dir/results/variants/candidateSV.vcf.gz > ${params.sample}.candidateSV.vcf
    gunzip -c ${params.sample}_manta_dir/results/variants/diploidSV.vcf.gz > ${params.sample}.diploidSV.vcf
    gunzip -c ${params.sample}_manta_dir/results/variants/candidateSmallIndels.vcf.gz > ${params.sample}.candidateSmallIndels.vcf
    rm -rf ${params.sample}_manta_dir/
    """

	
}
