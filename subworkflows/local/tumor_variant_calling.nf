//
// TUMOR VARIANT CALLING
// Should be only run on patients without normal sample
//

include { GATK_TUMOR_ONLY_SOMATIC_VARIANT_CALLING } from '../../subworkflows/nf-core/gatk_tumor_only_somatic_variantcalling/main'
include { FREEBAYES                               } from '../../modules/nf-core/modules/freebayes/main'
include { MANTA_TUMORONLY                         } from '../../modules/nf-core/modules/manta/tumor_only/main'
include { MSISENSOR_MSI                           } from '../../modules/nf-core/modules/msisensor/msi/main'
include { MSISENSOR_SCAN                          } from '../../modules/nf-core/modules/msisensor/scan/main'
include { STRELKA_GERMLINE as STRELKA_TUMORONLY   } from '../../modules/nf-core/modules/strelka/germline/main'

workflow TUMOR_ONLY_VARIANT_CALLING {
    take:
        tools             // Mandatory, list of tools to apply
        cram_recalibrated // channel: [mandatory] cram
        dbsnp             // channel: [mandatory] dbsnp
        dbsnp_tbi         // channel: [mandatory] dbsnp_tbi
        dict              // channel: [mandatory] dict
        fasta             // channel: [mandatory] fasta
        fasta_fai         // channel: [mandatory] fasta_fai
        intervals         // channel: [mandatory] intervals/target regions
        intervals_bed_gz_tbi // channel: [mandatory] intervals/target regions index zipped and indexed
        num_intervals     // val: number of intervals that are used to parallelize exection, either based on capture kit or GATK recommended for WGS
        //target_bed        // channel: [optional]  target_bed
        //target_bed_gz_tbi // channel: [optional]  target_bed_gz_tbi

    main:

    ch_versions = Channel.empty()

    cram_recalibrated.combine(intervals)
        .map{ meta, cram, crai, intervals ->
            new_meta = meta.clone()
            new_meta.id = meta.sample + "_" + intervals.simpleName
            [new_meta, cram, crai, intervals]
        }.set{cram_recalibrated_intervals}

    cram_recalibrated.combine(intervals_bed_gz_tbi)
        .map{ meta, cram, crai, bed, tbi ->
            new_meta = meta.clone()
            new_meta.id = meta.sample + "_" + bed.simpleName
            [new_meta, cram, crai, bed, tbi]
        }.set{cram_recalibrated_intervals_gz_tbi}

    if (tools.contains('mutect2')) {
    }

    if (tools.contains('freebayes')){
        //TODO: merge parallelized vcfs, index them all
        //TODO: research if multiple targets can be provided
        //TODO: Pass over dbsnp/knwon_indels?
        cram_recalibrated.combine(intervals)
        .map{ meta, cram, crai, intervals ->
            new_meta = meta.clone()
            new_meta.id = meta.sample + "_" + intervals.simpleName
            [new_meta, cram, crai, [], [], intervals]
        }.set{cram_recalibrated_intervals_freebayes}

        FREEBAYES(
            cram_recalibrated_intervals_freebayes,
            fasta,
            fasta_fai,
            [],
            [],
            []
        )
        freebayes_vcf_gz = FREEBAYES.out.vcf
        ch_versions = ch_versions.mix(FREEBAYES.out.versions)
    }

    if (tools.contains('manta')){
        //TODO: test data not running
        //TODO: merge parallelized vcfs, index them all
        //TODO: Pass over dbsnp/knwon_indels?
        cram_recalibrated
        .map{ meta, cram, crai ->
            new_meta = meta.clone()
            //new_meta.id = meta.sample + "_" + intervals.simpleName
            [new_meta, cram, crai, [], []]
        }.set{cram_recalibrated_manta}

        MANTA_TUMORONLY(
            cram_recalibrated_manta,
            fasta,
            fasta_fai
        )
        manta_candidate_small_indels_vcf_tbi = MANTA_TUMORONLY.out.candidate_small_indels_vcf.join(MANTA_TUMORONLY.out.candidate_small_indels_vcf_tbi)
        manta_candidate_sv_vcf_tbi = MANTA_TUMORONLY.out.candidate_sv_vcf.join(MANTA_TUMORONLY.out.candidate_sv_vcf_tbi)
        manta_diploid_sv_vcf_tbi = MANTA_TUMORONLY.out.diploid_sv_vcf.join(MANTA_TUMORONLY.out.diploid_sv_vcf)

        ch_versions = ch_versions.mix(MANTA_TUMORONLY.out.versions)
    }

    if (tools.contains('strelka')) {
        //TODO: research if multiple targets can be provided
        //TODO: merge parallelized vcfs, index them all
        //TODO: Pass over dbsnp/knwon_indels?

        STRELKA_TUMORONLY(
            cram_recalibrated_intervals_gz_tbi,
            fasta,
            fasta_fai
            )

        strelka_vcf_tbi = STRELKA_TUMORONLY.out.vcf.join(STRELKA_TUMORONLY.out.vcf_tbi)
        strelka_genome_vcf_tbi = STRELKA_TUMORONLY.out.genome_vcf.join(STRELKA_TUMORONLY.out.genome_vcf_tbi)

        ch_versions = ch_versions.mix(STRELKA_TUMORONLY.out.versions)

    }

    if (tools.contains('msisensor')){

    }

    emit:
    strelka_vcf_tbi        = Channel.empty()
    strelka_genome_vcf_tbi = Channel.empty()
}
