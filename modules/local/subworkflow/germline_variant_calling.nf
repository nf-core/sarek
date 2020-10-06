/*
================================================================================
                            GERMLINE VARIANT CALLING
================================================================================
*/

include { GATK_HAPLOTYPECALLER as HAPLOTYPECALLER } from '../../nf-core/software/gatk/haplotypecaller'
include { GATK_GENOTYPEGVCF as GENOTYPEGVCF }       from '../../nf-core/software/gatk/genotypegvcf'


workflow GERMLINE_VARIANT_CALLING {
    take:
        bam_variant_calling // channel: [mandatory] bam
        intervals           // channel: [mandatory] intervals
        tools               //   list:  [mandatory] list of tools
        target_bed          // channel: [optional]  target_bed
        dbsnp               // channel: [mandatory] dbsnp
        dbsnp_tbi           // channel: [mandatory] dbsnp_tbi
        fasta               // channel: [mandatory] fasta
        fai                 // channel: [mandatory] fai

    main:

    vcfGenotypeGVCFs = Channel.empty()

    if ('haplotypecaller' in tools) {
        bam_haplotypecaller = bam_variant_calling.combine(intervals)

        // STEP GATK HAPLOTYPECALLER.1

        HAPLOTYPECALLER(
            bam_haplotypecaller,
            dbsnp,
            dbsnp_tbi,
            dict,
            fasta,
            fai)

   
        // STEP GATK HAPLOTYPECALLER.2
        GENOTYPEGVCF(
            HAPLOTYPECALLER.out.gvcfGenotypeGVCFs,
            dbsnp,
            dbsnp_tbi,
            dict,
            fasta,
            fai)

        GENOTYPEGVCF.out.map{ name, meta, vcf -> 
            patient = meta.patient
            sample  = meta.sample
            gender  = meta.gender
            status  = meta.status
            [name, patient, sample, gender, status, vcf] 
        }.groupTuple(by: [0,1,2,])
         .set{ vcfGenotypeGVCFs }
    }

    emit:
        vcfGenotypeGVCFs = vcfGenotypeGVCFs
}
