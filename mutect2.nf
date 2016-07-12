#!/usr/bin/env nextflow

params.tumor_bam = "tcga.cl.tumor__1.recal.bam"
params.tumor_bam_idx = "tcga.cl.tumor__1.recal.bai"
params.normal_bam = "tcga.cl.normal__0.recal.bam"
params.normal_bam_idx = "tcga.cl.normal__0.recal.bai"

tumor_bam = file(params.tumor_bam)
tumor_bam_idx = file(params.tumor_bam_idx)
normal_bam = file(params.normal_bam)
normal_bam_idx = file(params.normal_bam_idx)

// define intervals file by --intervals
intervalsFile = file(params.intervals)
intervals = Channel
    .from(intervalsFile.readLines())

process runIntervals {
    input:
    file tumor_bam
    file tumor_bam_idx
    file normal_bam 
    file normal_bam_idx
    set genomicInterval from intervals

    output:
    file '*.vcf' into finalsVCFs

    """
    java -Xmx2g -jar /home/szilva/dev/mutect/GenomeAnalysisTK.jar \
    -T MuTect2 \
    -R /sw/data/uppnex/ToolBox/ReferenceAssemblies/hg38make/bundle/2.8/b37/human_g1k_v37_decoy.fasta \
    -I:tumor ${tumor_bam}\
    -I:normal ${normal_bam} \
    --dbsnp /sw/data/uppnex/ToolBox/ReferenceAssemblies/hg38make/bundle/2.8/b37/dbsnp_138.b37.vcf \
    --cosmic /sw/data/uppnex/ToolBox/ReferenceAssemblies/hg38make/bundle/2.8/b37/b37_cosmic_v54_120711.vcf \
    -L \"$genomicInterval\" \
    -o ${genomicInterval}.vcf
    """
}
