/*
================================================================================
                            GERMLINE VARIANT CALLING
================================================================================
*/

include { GATK_HAPLOTYPECALLER as HAPLOTYPECALLER } from '../../nf-core/software/gatk/haplotypecaller'
include { GATK_GENOTYPEGVCF as GENOTYPEGVCF }       from '../../nf-core/software/gatk/genotypegvcf'
include { CONCAT_VCF as CONCAT_GVCF;
          CONCAT_VCF as CONCAT_HAPLOTYPECALLER}     from '../process/concat_vcf'
include { STRELKA_GERMLINE }                        from '../../nf-core/software/strelka/germline'

workflow GERMLINE_VARIANT_CALLING {
    take:
        bam        // channel: [mandatory] bam
        intervals  // channel: [mandatory] intervals
        tools      //   list:  [mandatory] list of tools
        target_bed // channel: [optional]  target_bed
        dict       // channel: [mandatory] dict
        dbsnp      // channel: [mandatory] dbsnp
        dbsnp_tbi  // channel: [mandatory] dbsnp_tbi
        fasta      // channel: [mandatory] fasta
        fai        // channel: [mandatory] fai
        modules    //     map: [mandatory] maps for modules

    main:

    haplotypecaller_gvcf = Channel.empty()
    haplotypecaller_vcf  = Channel.empty()
    strelka_vcf          = Channel.empty()

    if ('haplotypecaller' in tools) {
        haplotypecaller_interval_bam = bam.combine(intervals)

        // STEP GATK HAPLOTYPECALLER.1

        HAPLOTYPECALLER(
            haplotypecaller_interval_bam,
            dbsnp,
            dbsnp_tbi,
            dict,
            fasta,
            fai)

        haplotypecaller_interval_gvcf = HAPLOTYPECALLER.out.gvcf.map{ meta, vcf ->
            patient = meta.patient
            sample  = meta.sample
            gender  = meta.gender
            status  = meta.status
            [ patient, sample, gender, status, vcf] 
        }.groupTuple(by: [0,1])

        haplotypecaller_interval_gvcf = haplotypecaller_interval_gvcf.map { patient, sample, gender, status, vcf ->
            def meta = [:]
            meta.patient = patient
            meta.sample  = sample
            meta.gender  = gender[0]
            meta.status  = status[0]
            meta.id      = meta.sample
            [ meta, vcf ]
        }

        CONCAT_GVCF(
            haplotypecaller_interval_gvcf,
            fai,
            target_bed,
            modules['concat_vcf_haplotypecaller_gvcf'])

        haplotypecaller_gvcf = CONCAT_GVCF.out.vcf

        // STEP GATK HAPLOTYPECALLER.2

        GENOTYPEGVCF(
            HAPLOTYPECALLER.out.interval_gvcf,
            dbsnp,
            dbsnp_tbi,
            dict,
            fasta,
            fai)

        haplotypecaller_interval_vcf = GENOTYPEGVCF.out.map{ meta, vcf ->
            patient = meta.patient
            sample  = meta.sample
            gender  = meta.gender
            status  = meta.status
            [ patient, sample, gender, status, vcf] 
        }.groupTuple(by: [0,1])

        haplotypecaller_interval_vcf = haplotypecaller_interval_vcf.map { patient, sample, gender, status, vcf ->
            def meta = [:]
            meta.patient = patient
            meta.sample  = sample
            meta.gender  = gender[0]
            meta.status  = status[0]
            meta.id      = meta.sample
            [ meta, vcf ]
        }

        CONCAT_HAPLOTYPECALLER(
            haplotypecaller_interval_vcf,
            fai,
            target_bed,
            modules['concat_vcf_haplotypecaller'])

        haplotypecaller_vcf = CONCAT_GVCF.out.vcf
    }

    if ('strelka' in tools) {
        STRELKA_GERMLINE(
            bam,
            fasta,
            fai,
            target_bed,
            modules['strelka_germline'])

        strelka_vcf = STRELKA_GERMLINE.out.vcf
    }

    emit:
        haplotypecaller_gvcf = haplotypecaller_gvcf
        haplotypecaller_vcf  = haplotypecaller_vcf
        strelka_vcf          = strelka_vcf
}
