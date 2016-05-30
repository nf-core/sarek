#!/usr/bin/env nextflow
//This script runs tumor normal analysis using delly2
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

params.delly2_path=""
params.SVDB_path=""

process delly_del{

    module 'bioinfo-tools'
    module 'bcftools'

    input:
    file genome_file

    file tumor_bam
    file normal_bam

    file tumor_bai
    file normal_bai

    output:
       file "${params.sample}.DELLY_del.vcf" into delly_del_vcf

    """
    ${params.delly2_path} call -t DEL -g {genome_file} ${normal_bam} ${tumor_bam} -o ${params.sample}.DELLY_del.bcf
    bcftools view ${params.sample}.DELLY_del.bcf > ${params.sample}.DELLY_del.vcf
    """

	
}

process delly_dup{

    module 'bioinfo-tools'
    module 'bcftools'


    input:
    file genome_file

    file tumor_bam
    file normal_bam

    file tumor_bai
    file normal_bai

    output:
       file "${params.sample}.DELLY_dup.vcf" into delly_dup_vcf

    """
    ${params.delly2_path} call -t DUP -g {genome_file} ${normal_bam} ${tumor_bam} -o ${params.sample}.DELLY_dup.bcf
    bcftools view ${params.sample}.DELLY_dup.bcf > ${params.sample}.DELLY_dup.vcf
    """
}
process delly_tra{

    module 'bioinfo-tools'
    module 'bcftools'

    input:
    file genome_file

    file tumor_bam
    file normal_bam

    file tumor_bai
    file normal_bai

    output:
       file "${params.sample}.DELLY_tra.vcf" into delly_tra_vcf

    """
    ${params.delly2_path} call -t TRA -g {genome_file} ${normal_bam} ${tumor_bam} -o ${params.sample}.DELLY_tra.bcf
    bcftools view ${params.sample}.DELLY_tra.bcf > ${params.sample}.DELLY_tra.vcf
    """
}
process delly_inv{

    module 'bioinfo-tools'
    module 'bcftools'

    input:
    file genome_file

    file tumor_bam
    file normal_bam

    file tumor_bai
    file normal_bai

    output:
       file "${params.sample}.DELLY_inv.vcf" into delly_inv_vcf

    """
    ${params.delly2_path} call -t INV -g {genome_file} ${normal_bam} ${tumor_bam} -o ${params.sample}.DELLY_inv.bcf
    bcftools view ${params.sample}.DELLY_inv.bcf > ${params.sample}.DELLY_inv.vcf
    """
}
process delly_ins{

    module 'bioinfo-tools'
    module 'bcftools'

    input:
    file genome_file

    file tumor_bam
    file normal_bam

    file tumor_bai
    file normal_bai

    output:
       file "${params.sample}.DELLY_ins.vcf" into delly_ins_vcf

    """
    ${params.delly2_path} call -t INS -g {genome_file} ${normal_bam} ${tumor_bam} -o ${params.sample}.DELLY_ins.bcf
    bcftools view ${params.sample}.DELLY_ins.bcf > ${params.sample}.DELLY_ins.vcf
    """
}

//process delly_merge{
//    input:
//
//    file delly_inv_vcf
//    file delly_ins_vcf
//    file delly_tra_vcf
 //   file delly_dup_vcf
//    file delly_del_vcf
//
//    file delly_tumor_inv_vcf
//    file delly_tumor_ins_vcf
//    file delly_tumor_tra_vcf
//    file delly_tumor_dup_vcf
//    file delly_tumor_del_vcf

//    output:
//        file "${params.sample}.DELLY_tumor.vcf" into delly_tumor_vcf
//        file "${params.sample}.DELLY_normal.vcf" into delly_normal_vcf




//}


