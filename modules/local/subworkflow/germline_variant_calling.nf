/*
================================================================================
                            GERMLINE VARIANT CALLING
================================================================================
*/

include { GATK_HAPLOTYPECALLER as HAPLOTYPECALLER } from '../../nf-core/software/gatk/haplotypecaller'
include { GATK_GENOTYPEGVCF as GENOTYPEGVCF }       from '../../nf-core/software/gatk/genotypegvcf'
include { CONCAT_VCF as CONCAT_GVCF;
          CONCAT_VCF as CONCAT_HAPLOTYPECALLER}     from '../process/concat_vcf'


workflow GERMLINE_VARIANT_CALLING {
    take:
        bam_variant_calling // channel: [mandatory] bam
        intervals           // channel: [mandatory] intervals
        tools               //   list:  [mandatory] list of tools
        target_bed          // channel: [optional]  target_bed
        dict                // channel: [mandatory] dict
        dbsnp               // channel: [mandatory] dbsnp
        dbsnp_tbi           // channel: [mandatory] dbsnp_tbi
        fasta               // channel: [mandatory] fasta
        fai                 // channel: [mandatory] fai
        modules             //     map: [mandatory] maps for modules

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

        haplotypecallergvcf = HAPLOTYPECALLER.out.gvcf.map{ meta, vcf ->
            patient = meta.patient
            sample  = meta.sample
            gender  = meta.gender
            status  = meta.status
            [ patient, sample, gender, status, vcf] 
        }.groupTuple(by: [0,1])

        haplotypecallergvcf = haplotypecallergvcf.map { patient, sample, gender, status, vcf ->
            def meta = [:]
            meta.patient = patient
            meta.sample  = sample
            meta.gender  = gender[0]
            meta.status  = status[0]
            meta.id      = meta.sample
            [ meta, vcf ]
        }

        CONCAT_GVCF(
            haplotypecallergvcf,
            fai,
            target_bed,
            modules['concat_vcf_haplotypecallergvcf'])
   
        // STEP GATK HAPLOTYPECALLER.2

        GENOTYPEGVCF(
            HAPLOTYPECALLER.out.interval_gvcf,
            dbsnp,
            dbsnp_tbi,
            dict,
            fasta,
            fai)

        haplotypecallervcf = GENOTYPEGVCF.out.map{ meta, vcf ->
            patient = meta.patient
            sample  = meta.sample
            gender  = meta.gender
            status  = meta.status
            [ patient, sample, gender, status, vcf] 
        }.groupTuple(by: [0,1])

        haplotypecallervcf = haplotypecallervcf.map { patient, sample, gender, status, vcf ->
            def meta = [:]
            meta.patient = patient
            meta.sample  = sample
            meta.gender  = gender[0]
            meta.status  = status[0]
            meta.id      = meta.sample
            [ meta, vcf ]
        }

        CONCAT_HAPLOTYPECALLER(
            haplotypecallervcf,
            fai,
            target_bed,
            modules['concat_vcf_haplotypecaller'])

    }

    emit:
        vcfGenotypeGVCFs = vcfGenotypeGVCFs
}
