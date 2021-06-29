/*
========================================================================================
    GERMLINE VARIANT CALLING
========================================================================================
*/

params.haplotypecaller_options        = [:]
params.genotypegvcf_options           = [:]
params.concat_gvcf_options            = [:]
params.concat_haplotypecaller_options = [:]
params.strelka_options                = [:]

include { GATK4_HAPLOTYPECALLER as HAPLOTYPECALLER } from '../../modules/nf-core/software/gatk4/haplotypecaller/main' addParams(options: params.haplotypecaller_options)
include { GATK4_GENOTYPEGVCF as GENOTYPEGVCF }       from '../../modules/nf-core/software/gatk4/genotypegvcf/main'    addParams(options: params.genotypegvcf_options)
include { CONCAT_VCF as CONCAT_GVCF }                from '../../modules/local/concat_vcf/main'                       addParams(options: params.concat_gvcf_options)
include { CONCAT_VCF as CONCAT_HAPLOTYPECALLER }     from '../../modules/local/concat_vcf/main'                       addParams(options: params.concat_haplotypecaller_options)
include { STRELKA_GERMLINE as STRELKA }              from '../../modules/nf-core/software/strelka/germline/main'      addParams(options: params.strelka_options)

workflow GERMLINE_VARIANT_CALLING {
    take:
        tools
        cram              // channel: [mandatory] cram
        dbsnp             // channel: [mandatory] dbsnp
        dbsnp_tbi         // channel: [mandatory] dbsnp_tbi
        dict              // channel: [mandatory] dict
        fai               // channel: [mandatory] fai
        fasta             // channel: [mandatory] fasta
        intervals         // channel: [mandatory] intervals
        num_intervals
        target_bed        // channel: [optional]  target_bed
        target_bed_gz_tbi // channel: [optional]  target_bed_gz_tbi

    main:

    haplotypecaller_gvcf = Channel.empty()
    haplotypecaller_vcf  = Channel.empty()
    strelka_vcf          = Channel.empty()

    no_intervals = false
    if (intervals == []) no_intervals = true

    if ('haplotypecaller' in tools) {

        cram.combine(intervals).map{ meta, cram, crai, intervals ->
            new_meta = meta.clone()
            new_meta.id = meta.sample + "_" + intervals.baseName
            [new_meta, cram, crai, intervals]
        }.set{haplotypecaller_interval_cram}

        // STEP GATK HAPLOTYPECALLER.1
        HAPLOTYPECALLER(
            haplotypecaller_interval_cram,
            dbsnp,
            dbsnp_tbi,
            dict,
            fasta,
            fai,
            no_intervals)

        haplotypecaller_raw = HAPLOTYPECALLER.out.vcf.map{ meta,vcf ->
            meta.id = meta.sample
            [meta, vcf]
        }.groupTuple(size: num_intervals)

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

    if ('strelka' in tools) {
        cram.dump()
        fasta.dump()
        fai.dump()
        target_bed_gz_tbi.dump(tag:'target_bed_gz_tbi')
        STRELKA(
            cram,
            fasta,
            fai,
            target_bed_gz_tbi)

        strelka_vcf = STRELKA.out.vcf
    }

    //TODO add all the remaining variant caller:
    //FREEBAYES
    //SAMTOOLS MPILEUP
    //TIDDIT
    //MANTA

    emit:
        haplotypecaller_gvcf = haplotypecaller_gvcf
        haplotypecaller_vcf  = haplotypecaller_vcf
        strelka_vcf          = strelka_vcf
}
