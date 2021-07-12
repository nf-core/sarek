/*
========================================================================================
    GERMLINE VARIANT CALLING
========================================================================================
*/

params.haplotypecaller_options        = [:]
params.deepvariant_options            = [:]
params.genotypegvcf_options           = [:]
params.concat_gvcf_options            = [:]
params.concat_haplotypecaller_options = [:]
params.strelka_options                = [:]

include { GATK4_HAPLOTYPECALLER as HAPLOTYPECALLER } from '../../modules/nf-core/software/gatk4/haplotypecaller/main' addParams(options: params.haplotypecaller_options)
include { DEEPVARIANT }                              from '../../modules/local/deepvariant/main'                      addParams(options: params.deepvariant_options)
include { GATK4_GENOTYPEGVCF as GENOTYPEGVCF }       from '../../modules/nf-core/software/gatk4/genotypegvcf/main'    addParams(options: params.genotypegvcf_options)
include { CONCAT_VCF as CONCAT_GVCF }                from '../../modules/local/concat_vcf/main'                       addParams(options: params.concat_gvcf_options)
include { CONCAT_VCF as CONCAT_HAPLOTYPECALLER }     from '../../modules/local/concat_vcf/main'                       addParams(options: params.concat_haplotypecaller_options)
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

    main:

    haplotypecaller_gvcf = Channel.empty()
    haplotypecaller_vcf  = Channel.empty()
    strelka_vcf          = Channel.empty()
    no_intervals = false

    if (intervals == []) no_intervals = true

    if ('haplotypecaller' in params.tools.toLowerCase()) {

        haplotypecaller_interval_bam = bam.combine(intervals)

        bam.combine(intervals).map{ meta, bam, bai, intervals ->
            new_meta = meta.clone()
            new_meta.id = meta.sample + "_" + intervals.baseName
            [new_meta, bam, bai, intervals]
        }.set{haplotypecaller_interval_bam}

        // STEP GATK HAPLOTYPECALLER.1

        HAPLOTYPECALLER(
            haplotypecaller_interval_bam,
            dbsnp,
            dbsnp_tbi,
            dict,
            fasta,
            fai,
            no_intervals)

        haplotypecaller_raw = HAPLOTYPECALLER.out.vcf.map{ meta,vcf ->
            meta.id = meta.sample
            [meta, vcf]
        }.groupTuple()

        CONCAT_GVCF(
            haplotypecaller_raw,
            fai,
            target_bed)

        haplotypecaller_gvcf = CONCAT_GVCF.out.vcf

        // STEP GATK HAPLOTYPECALLER.2

        GENOTYPEGVCF(
            HAPLOTYPECALLER.out.interval_vcf,
            dbsnp,
            dbsnp_tbi,
            dict,
            fasta,
            fai,
            no_intervals)

        haplotypecaller_results = GENOTYPEGVCF.out.vcf.map{ meta, vcf ->
            meta.id = meta.sample
            [meta, vcf]
        }.groupTuple()

        CONCAT_HAPLOTYPECALLER(
            haplotypecaller_results,
            fai,
            target_bed)
        
        haplotypecaller_vcf = CONCAT_HAPLOTYPECALLER.out.vcf
    }

    if ('deepvariant' in params.tools.toLowerCase()) {

        DEEPVARIANT(
            bam,
            fasta,
            fai)

        haplotypecaller_vcf = DEEPVARIANT.out.vcf
        haplotypecaller_gvcf = DEEPVARIANT.out.gvcf
    }

    if ('strelka' in params.tools.toLowerCase()) {
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
