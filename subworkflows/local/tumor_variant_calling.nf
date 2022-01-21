//
// TUMOR VARIANT CALLING
// Should be only run on patients without normal sample
//


include { BGZIP as BGZIP_FREEBAYES                    } from '../../modules/local/bgzip'
include { BGZIP as BGZIP_MANTA_SMALL_INDELS           } from '../../modules/local/bgzip'
include { BGZIP as BGZIP_MANTA_SV                     } from '../../modules/local/bgzip'
include { BGZIP as BGZIP_MANTA_DIPLOID                } from '../../modules/local/bgzip'
include { BGZIP as BGZIP_STRELKA                      } from '../../modules/local/bgzip'
include { CONCAT_VCF as CONCAT_VCF_FREEBAYES          } from '../../modules/local/concat_vcf/main'
include { CONCAT_VCF as CONCAT_VCF_MANTA_SMALL_INDELS } from '../../modules/local/concat_vcf/main'
include { CONCAT_VCF as CONCAT_VCF_MANTA_SV           } from '../../modules/local/concat_vcf/main'
include { CONCAT_VCF as CONCAT_VCF_MANTA_DIPLOID      } from '../../modules/local/concat_vcf/main'
include { CONCAT_VCF as CONCAT_VCF_STRELKA            } from '../../modules/local/concat_vcf/main'
include { FREEBAYES                                   } from '../../modules/nf-core/modules/freebayes/main'
include { GATK_TUMOR_ONLY_SOMATIC_VARIANT_CALLING     } from '../../subworkflows/nf-core/gatk_tumor_only_somatic_variant_calling/main'
include { MANTA_TUMORONLY                             } from '../../modules/nf-core/modules/manta/tumoronly/main'
include { STRELKA_GERMLINE as STRELKA_TUMORONLY       } from '../../modules/nf-core/modules/strelka/germline/main'

workflow TUMOR_ONLY_VARIANT_CALLING {
    take:
        tools                   // Mandatory, list of tools to apply
        cram_recalibrated       // channel: [mandatory] cram
        dbsnp                   // channel: [mandatory] dbsnp
        dbsnp_tbi               // channel: [mandatory] dbsnp_tbi
        dict                    // channel: [mandatory] dict
        fasta                   // channel: [mandatory] fasta
        fasta_fai               // channel: [mandatory] fasta_fai
        intervals               // channel: [mandatory] intervals/target regions
        intervals_bed_gz_tbi    // channel: [mandatory] intervals/target regions index zipped and indexed
        num_intervals           // val: number of intervals that are used to parallelize exection, either based on capture kit or GATK recommended for WGS
        germline_resource
        germline_resource_tbi // channel
        panel_of_normals
        panel_of_normals_tbi
        //target_bed        // channel: [optional]  target_bed
        //target_bed_gz_tbi // channel: [optional]  target_bed_gz_tbi

    main:

    ch_versions                 = Channel.empty()
    freebayes_vcf_gz_tbi        = Channel.empty()
    strelka_vcf_gz_tbi          = Channel.empty()

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

    if (tools.contains('freebayes')){
        //TODO: Pass over dbsnp/knwon_indels?
        cram_recalibrated.combine(intervals).map{ meta, cram, crai, intervals ->
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

        BGZIP_FREEBAYES(FREEBAYES.out.vcf)
        freebayes_vcf_to_concat = BGZIP_FREEBAYES.out.vcf.groupTuple(size: num_intervals)

        CONCAT_VCF_FREEBAYES(freebayes_vcf_to_concat,fasta_fai, intervals)
        freebayes_vcf_gz_tbi = CONCAT_VCF_FREEBAYES.out.vcf

        ch_versions = ch_versions.mix(FREEBAYES.out.versions)
        ch_versions = ch_versions.mix(BGZIP_FREEBAYES.out.versions)
        ch_versions = ch_versions.mix(CONCAT_VCF_FREEBAYES.out.versions)
    }

    if (tools.contains('mutect2')) {
        // GATK_TUMOR_ONLY_SOMATIC_VARIANT_CALLING(
        //     cram_recalibrated,
        //     fasta,
        //     fasta_fai,
        //     dict,
        //     germline_resource,
        //     germline_resource_tbi,
        //     panel_of_normals,
        //     panel_of_normals_tbi,
        //     intervals
        // )
        //TODO: mutectSTATS
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
        manta_tumor_sv_vcf_tbi = MANTA_TUMORONLY.out.tumor_sv_vcf.join(MANTA_TUMORONLY.out.tumor_sv_vcf)

        ch_versions = ch_versions.mix(MANTA_TUMORONLY.out.versions)
    }

    if (tools.contains('strelka')) {
        //TODO: research if multiple targets can be provided: waiting for reply
        //TODO: Pass over dbsnp/knwon_indels?

        STRELKA_TUMORONLY(
            cram_recalibrated_intervals_gz_tbi,
            fasta,
            fasta_fai
            )

        BGZIP_STRELKA(STRELKA_TUMORONLY.out.vcf)
        strelka_vcf_to_concat = BGZIP_STRELKA.out.vcf.groupTuple(size: num_intervals)

        CONCAT_VCF_STRELKA(strelka_vcf_to_concat,fasta_fai, intervals)
        strelka_vcf_gz_tbi = CONCAT_VCF_STRELKA.out.vcf

        ch_versions = ch_versions.mix(STRELKA_TUMORONLY.out.versions)
        ch_versions = ch_versions.mix(BGZIP_STRELKA.out.versions)
        ch_versions = ch_versions.mix(CONCAT_VCF_STRELKA.out.versions)
    }

    //if (tools.contains('msisensor')){
        //TODO: Requires a baseline optimally per tumor/sequencing type. I am not convinced it is in
        // our scope to tupport this at this time
        //MSISENSORPRO_MSI(cram_recalibrated, msisensorpro_scan)
    //}

    emit:
    strelka_vcf_tbi        = Channel.empty()
    strelka_genome_vcf_tbi = Channel.empty()
}
