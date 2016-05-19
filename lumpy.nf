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


process normal_discordant{

    module 'bioinfo-tools'
    module 'samtools'

    input:
    
    file normal_bam
    file normal_bai

    output:
    file "${params.sample}_normal.D.bam" into normal_discordant_bam

    """
    samtools view -b -F 1294 ${params.normal_bam} > ${params.sample}_normal.D.unsorted.bam
    samtools sort ${params.sample}_normal.D.unsorted.bam ${params.sample}_normal.D
    rm ${params.sample}_normal.D.unsorted.bam
    """

}

process tumor_discordant{

    module 'bioinfo-tools'
    module 'samtools'

    input:
    file tumor_bam
    file tumor_bai

    output:
    file "${params.sample}_tumor.D.bam" into tumor_discordant_bam

    """
    samtools view -b -F 1294 ${params.tumor_bam} > ${params.sample}_tumor.D.unsorted.bam
    samtools sort ${params.sample}_tumor.D.unsorted.bam ${params.sample}_tumor.D
    rm ${params.sample}_tumor.D.unsorted.bam 
    """

}

process normal_split{

    module 'bioinfo-tools'
    module 'samtools'

    input:
    file normal_bam
    file normal_bai

    output:
    file "${params.sample}_normal.S.bam" into normal_split_bam

    """
    samtools view -h ${params.normal_bam} | ${params.splitreads_extract} -i stdin | samtools view -Sb - > ${params.sample}_normal.S.unsorted.bam
    samtools sort ${params.sample}_normal.S.unsorted.bam ${params.sample}_normal.S
    rm ${params.sample}_normal.S.unsorted.bam
    """
}

process tumor_split{

    module 'bioinfo-tools'
    module 'samtools'

    input:
    file tumor_bam
    file tumor_bai

    output:
    file "${params.sample}_tumor.S.bam" into tumor_split_bam


    """
    samtools view -h ${params.tumor_bam} | ${params.splitreads_extract} -i stdin | samtools view -Sb - > ${params.sample}_tumor.S.unsorted.bam
    samtools sort ${params.sample}_tumor.S.unsorted.bam ${params.sample}_tumor.S
    rm ${params.sample}_tumor.S.unsorted.bam
    """

}

process lumpy{

    module 'bioinfo-tools'
    module 'LUMPY/0.2.12'

    input:

    file tumor_bam
    file normal_bam

    file tumor_bai
    file normal_bai

    file normal_discordant_bam
    file tumor_discordant_bam

    file normal_split_bam
    file tumor_split_bam

    output:
       file "${params.sample}.tumor_normal.lumpySV.vcf" into lumpy_vcf

    """
    lumpyexpress -B ${params.normal_bam},${params.tumor_bam} -S ${normal_split_bam},${tumor_split_bam} -D ${normal_discordant_bam},${tumor_discordant_bam} -o ${params.sample}.tumor_normal.lumpySV.vcf
    rm ${normal_split_bam} ${tumor_split_bam} ${normal_discordant_bam} ${tumor_discordant_bam}
    """
}
