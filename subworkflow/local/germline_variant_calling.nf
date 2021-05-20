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

include { GATK4_HAPLOTYPECALLER as HAPLOTYPECALLER } from '../../modules/nf-core/software/gatk4/haplotypecaller/main' addParams(options: params.haplotypecaller_options)
include { GATK4_GENOTYPEGVCF as GENOTYPEGVCF }       from '../../modules/nf-core/software/gatk4/genotypegvcf/main'    addParams(options: params.genotypegvcf_options)
include { CONCAT_VCF as CONCAT_GVCF }                from '../../modules/local/concat_vcf'                            addParams(options: params.concat_gvcf_options)
include { CONCAT_VCF as CONCAT_HAPLOTYPECALLER }     from '../../modules/local/concat_vcf'                            addParams(options: params.concat_haplotypecaller_options)
include { STRELKA_GERMLINE as STRELKA }              from '../../modules/nf-core/software/strelka/germline/main'      addParams(options: params.strelka_options)

workflow GERMLINE_VARIANT_CALLING {
    take:
        bam               // channel: [mandatory] bam
        dbsnp             // channel: [mandatory] dbsnp
        dbsnp_tbi         // channel: [mandatory] dbsnp_tbi
        dict              // channel: [mandatory] dict
        fai               // channel: [mandatory] fai
        fasta             // channel: [mandatory] fasta
        intervals         // channel: [mandatory] intervals
        target_bed        // channel: [optional]  target_bed
        target_bed_gz_tbi // channel: [optional]  target_bed_gz_tbi
        tools             //   list:  [mandatory] list of tools

    main:

    haplotypecaller_gvcf = Channel.empty()
    haplotypecaller_vcf  = Channel.empty()
    strelka_vcf          = Channel.empty()
    no_intervals = false

    if (intervals == []) no_intervals = true

    if ('haplotypecaller' in tools) {
        haplotypecaller_interval_bam = bam.combine(intervals)

        // STEP GATK HAPLOTYPECALLER.1

        HAPLOTYPECALLER(
            haplotypecaller_interval_bam,
            dbsnp,
            dbsnp_tbi,
            dict,
            fasta,
            fai,
            no_intervals)

        haplotypecaller_gvcf = HAPLOTYPECALLER.out.gvcf.map{ meta, vcf ->
            patient = meta.patient
            sample  = meta.sample
            gender  = meta.gender
            status  = meta.status
            [ patient, sample, gender, status, vcf] 
        }.groupTuple(by: [0,1])

        haplotypecaller_gvcf = haplotypecaller_gvcf.map { patient, sample, gender, status, vcf ->
            def meta = [:]
            meta.patient = patient
            meta.sample  = sample
            meta.gender  = gender[0]
            meta.status  = status[0]
            meta.id      = meta.sample
            [ meta, vcf ]
        }

        CONCAT_GVCF(
            haplotypecaller_gvcf,
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
            fai,
            no_intervals)

        haplotypecaller_interval_vcf = GENOTYPEGVCF.out.vcf.map{ meta, vcf ->
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
            target_bed_gz_tbi)

        strelka_vcf = STRELKA.out.vcf
    }

    emit:
        haplotypecaller_gvcf = haplotypecaller_gvcf
        haplotypecaller_vcf  = haplotypecaller_vcf
        strelka_vcf          = strelka_vcf
}
