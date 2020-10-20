/*
================================================================================
                            GERMLINE VARIANT CALLING
================================================================================
*/

params.haplotypecaller_options        = [:]
params.genotypegvcf_options           = [:]
params.concat_gvcf_options            = [:]
params.concat_haplotypecaller_options = [:]
params.strelka_options                = [:]

include { GATK_HAPLOTYPECALLER as HAPLOTYPECALLER } from '../../nf-core/software/gatk/haplotypecaller' addParams(options: params.haplotypecaller_options)
include { GATK_GENOTYPEGVCF as GENOTYPEGVCF }       from '../../nf-core/software/gatk/genotypegvcf'    addParams(options: params.genotypegvcf_options)
include { CONCAT_VCF as CONCAT_GVCF }               from '../process/concat_vcf'                       addParams(options: params.concat_gvcf_options)
include { CONCAT_VCF as CONCAT_HAPLOTYPECALLER }    from '../process/concat_vcf'                       addParams(options: params.concat_haplotypecaller_options)
include { STRELKA_GERMLINE as STRELKA }             from '../../nf-core/software/strelka/germline'     addParams(options: params.strelka_options)

workflow GERMLINE_VARIANT_CALLING {
    take:
        bam        // channel: [mandatory] bam
        dbsnp      // channel: [mandatory] dbsnp
        dbsnp_tbi  // channel: [mandatory] dbsnp_tbi
        dict       // channel: [mandatory] dict
        fai        // channel: [mandatory] fai
        fasta      // channel: [mandatory] fasta
        intervals  // channel: [mandatory] intervals
        target_bed // channel: [optional]  target_bed
        tools      //   list:  [mandatory] list of tools

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
            target_bed)

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
            target_bed)

        haplotypecaller_vcf = CONCAT_GVCF.out.vcf
    }

    if ('strelka' in tools) {
        STRELKA(
            bam,
            fasta,
            fai,
            target_bed)

        strelka_vcf = STRELKA.out.vcf
    }

    emit:
        haplotypecaller_gvcf = haplotypecaller_gvcf
        haplotypecaller_vcf  = haplotypecaller_vcf
        strelka_vcf          = strelka_vcf
}
