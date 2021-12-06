//
// GERMLINE VARIANT CALLING
//

params.haplotypecaller_options        = [:]
params.deepvariant_options            = [:]
params.genotypegvcf_options           = [:]
params.concat_gvcf_options            = [:]
params.concat_haplotypecaller_options = [:]
params.strelka_options                = [:]

include { GATK4_HAPLOTYPECALLER as HAPLOTYPECALLER } from '../../modules/local/gatk4/haplotypecaller/main'      addParams(options: params.haplotypecaller_options)
include { DEEPVARIANT }                              from '../../modules/local/deepvariant/main'                addParams(options: params.deepvariant_options)
include { GATK4_GENOTYPEGVCF as GENOTYPEGVCF }       from '../../modules/local/gatk4/genotypegvcf/main'         addParams(options: params.genotypegvcf_options)
include { CONCAT_VCF as CONCAT_GVCF }                from '../../modules/local/concat_vcf/main'                 addParams(options: params.concat_gvcf_options)
include { CONCAT_VCF as CONCAT_HAPLOTYPECALLER }     from '../../modules/local/concat_vcf/main'                 addParams(options: params.concat_haplotypecaller_options)
include { STRELKA_GERMLINE as STRELKA }              from '../../modules/nf-core/modules/strelka/germline/main' addParams(options: params.strelka_options)

workflow GERMLINE_VARIANT_CALLING {
    take:
        tools
        cram              // channel: [mandatory] cram
        dbsnp             // channel: [mandatory] dbsnp
        dbsnp_tbi         // channel: [mandatory] dbsnp_tbi
        dict              // channel: [mandatory] dict
        fasta             // channel: [mandatory] fasta
        fasta_fai         // channel: [mandatory] fasta_fai
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
            fasta_fai)

        haplotypecaller_raw = HAPLOTYPECALLER.out.vcf.map{ meta,vcf ->
            meta.id = meta.sample
            [meta, vcf]
        }.groupTuple(size: num_intervals)

        CONCAT_GVCF(
            haplotypecaller_raw,
            fasta_fai,
            target_bed)

        haplotypecaller_gvcf = CONCAT_GVCF.out.vcf

        // STEP GATK HAPLOTYPECALLER.2

        HAPLOTYPECALLER.out.view()
        GENOTYPEGVCF(
            HAPLOTYPECALLER.out.interval_vcf,
            dbsnp,
            dbsnp_tbi,
            dict,
            fasta,
            fasta_fai)

        haplotypecaller_results = GENOTYPEGVCF.out.vcf.map{ meta, vcf ->
            meta.id = meta.sample
            [meta, vcf]
        }.groupTuple()

        CONCAT_HAPLOTYPECALLER(
            haplotypecaller_results,
            fasta_fai,
            target_bed)

        haplotypecaller_vcf = CONCAT_HAPLOTYPECALLER.out.vcf
    }

    if ('deepvariant' in tools) {

        DEEPVARIANT(
            bam,
            fasta,
            fasta_fai)

        deepvariant_vcf = DEEPVARIANT.out.vcf
        deepvariant_gvcf = DEEPVARIANT.out.gvcf
    }

    if ('strelka' in tools) {
        STRELKA(
            cram,
            fasta,
            fasta_fai,
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
