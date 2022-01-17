//
// GERMLINE VARIANT CALLING
//

include { GATK_JOINT_GERMLINE_VARIANT_CALLING } from '../../subworkflows/nf-core/joint_germline_variant_calling/main'
include { DEEPVARIANT                         } from '../../modules/nf-core/modules/deepvariant/main'
include { STRELKA_GERMLINE as STRELKA         } from '../../modules/nf-core/modules/strelka/germline/main'
include { FREEBAYES                           } from '../../modules/nf-core/modules/freebayes/main'
include { MANTA_GERMLINE                      } from '../../modules/nf-core/modules/manta/germline/main'
include { TIDDIT_SV                           } from '../../modules/nf-core/modules/tiddit/sv/main'
include { SAMTOOLS_MPILEUP                    } from '../../modules/nf-core/modules/samtools/mpileup/main'

workflow GERMLINE_VARIANT_CALLING {
    take:
        tools             // Mandatory, list of tools to apply
        cram              // channel: [mandatory] cram
        dbsnp             // channel: [mandatory] dbsnp
        dbsnp_tbi         // channel: [mandatory] dbsnp_tbi
        dict              // channel: [mandatory] dict
        fasta             // channel: [mandatory] fasta
        fasta_fai         // channel: [mandatory] fasta_fai
        intervals         // channel: [mandatory] intervals
        num_intervals     // val: number of intervals that are used to parallelize exection, either based on capture kit or GATK recommended for WGS
        target_bed        // channel: [optional]  target_bed
        target_bed_gz_tbi // channel: [optional]  target_bed_gz_tbi

    main:

    // haplotypecaller_gvcf = Channel.empty()
    // haplotypecaller_vcf  = Channel.empty()
    // strelka_vcf          = Channel.empty()

    // no_intervals = false
    // if (intervals == []) no_intervals = true

    if ('haplotypecaller' in tools) {
        //GATK_JOINT_GERMLINE_VARIANT_CALLING(

        //)

    //     cram.combine(intervals).map{ meta, cram, crai, intervals ->
    //         new_meta = meta.clone()
    //         new_meta.id = meta.sample + "_" + intervals.baseName
    //         [new_meta, cram, crai, intervals]
    //     }.set{haplotypecaller_interval_cram}

    //     // STEP GATK HAPLOTYPECALLER.1
    //     HAPLOTYPECALLER(
    //         haplotypecaller_interval_cram,
    //         dbsnp,
    //         dbsnp_tbi,
    //         dict,
    //         fasta,
    //         fasta_fai,
    //         no_intervals)

    //     haplotypecaller_raw = HAPLOTYPECALLER.out.vcf.map{ meta,vcf ->
    //         meta.id = meta.sample
    //         [meta, vcf]
    //     }.groupTuple(size: num_intervals)

    //     CONCAT_GVCF(
    //         haplotypecaller_raw,
    //         fasta_fai,
    //         target_bed)

    //     haplotypecaller_gvcf = CONCAT_GVCF.out.vcf

    //     // STEP GATK HAPLOTYPECALLER.2

    //     GENOTYPEGVCF(
    //         HAPLOTYPECALLER.out.interval_vcf,
    //         dbsnp,
    //         dbsnp_tbi,
    //         dict,
    //         fasta,
    //         fasta_fai,
    //         no_intervals)

    //     haplotypecaller_results = GENOTYPEGVCF.out.vcf.map{ meta, vcf ->
    //         meta.id = meta.sample
    //         [meta, vcf]
    //     }.groupTuple()

    //     CONCAT_HAPLOTYPECALLER(
    //         haplotypecaller_results,
    //         fasta_fai,
    //         target_bed)

    //     haplotypecaller_vcf = CONCAT_HAPLOTYPECALLER.out.vcf
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
            target_bed,
            target_bed_gz_tbi)

        strelka_vcf = STRELKA.out.vcf
    }

    if ('freebayes' in tools){
        FREEBAYES(
            cram,
            fasta,
            fasta_fai,
            target_bed,
            [],
            [],
            []
        )
    }

    if ('mpileup' in tools){
        SAMTOOLS_MPILEUP(
            cram,
            fasta
        )
    }

    if ('tiddit' in tools){
        TIDDIT_SV(
            cram,
            fasta,
            fasta_fai
        )
    }

    if ('manta' in tools){
        MANTA_GERMLINE(
            cram,
            fasta,
            fasta_fai,
            target_bed,
            target_bed_gz_tbi
        )
    }

    emit:
        haplotypecaller_gvcf = haplotypecaller_gvcf
        haplotypecaller_vcf  = haplotypecaller_vcf
        strelka_vcf          = strelka_vcf
}
