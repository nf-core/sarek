//
// PAIRED VARIANT CALLING
//
include { GATK_TUMOR_NORMAL_SOMATIC_VARIANT_CALLING     } from '../../subworkflows/nf-core/gatk_tumor_normal_somatic_variant_calling/main'
//include { GATK4_MUTECT2_SOMATIC as MUTECT2              } from '../../modules/nf-core/software/gatk4/mutect2/somatic/main'
include { MANTA_SOMATIC                                 } from '../../modules/nf-core/modules/manta/somatic/main'
include { MSISENSORPRO_MSI_SOMATIC                              } from '../../modules/nf-core/modules/msisensorpro/msi_somatic/main'
include { STRELKA_SOMATIC                               } from '../../modules/nf-core/modules/strelka/somatic/main'
include { STRELKA_SOMATIC as STRELKA_BP                 } from '../../modules/nf-core/modules/strelka/somatic/main'

workflow PAIR_VARIANT_CALLING {
    take:
        tools
        cram                  // channel: [mandatory] cram
        dbsnp                 // channel: [mandatory] dbsnp
        dbsnp_tbi             // channel: [mandatory] dbsnp_tbi
        dict                  // channel: [mandatory] dict
        fasta                 // channel: [mandatory] fasta
        fasta_fai             // channel: [mandatory] fasta_fai
        intervals               // channel: [mandatory] intervals/target regions
        intervals_bed_gz_tbi    // channel: [mandatory] intervals/target regions index zipped and indexed
        intervals_bed_combined_gz_tbi    // channel: [mandatory] intervals/target regions index zipped and indexed
        num_intervals           // val: number of intervals that are used to parallelize exection, either based on capture kit or GATK recommended for WGS
        msisensorpro_scan     // channel: [optional]  msisensorpro_scan
        // target_bed            // channel: [optional]  target_bed
        // target_bed_gz_tbi     // channel: [optional]  target_bed_gz_tbi
        germline_resource     // channel: [optional]  germline_resource
        germline_resource_tbi // channel: [optional]  germline_resource_tbi
        panel_of_normals      // channel: [optional]  panel_of_normals
        panel_of_normals_tbi  // channel: [optional]  panel_of_normals_tbi

    main:

    // cram.map{ meta, cram, crai ->
    //     patient = meta.patient
    //     sample  = meta.sample
    //     gender  = meta.gender
    //     status  = meta.status
    //     [patient, sample, gender, status, cram, crai]
    // }.branch{
    //     normal: it[3] == 0
    //     tumor:  it[3] == 1
    // }.set{ cram_to_cross }

    // cram_pair = cram_to_cross.normal.cross(cram_to_cross.tumor).map { normal, tumor ->
    //     def meta = [:]
    //     meta.patient = normal[0]
    //     meta.normal  = normal[1]
    //     meta.tumor   = tumor[1]
    //     meta.gender  = normal[2]
    //     meta.id      = "${meta.tumor}_vs_${meta.normal}".toString()

    //     [meta, normal[4], normal[5], tumor[4], tumor[5]]
    // }

    cram.combine(intervals).map{ meta, normal_cram, normal_crai, tumor_cram, tumor_crai, intervals ->
            new_meta = meta.clone()
            new_meta.id = meta.id + "_" + intervals.baseName
            [new_meta, normal_cram, normal_crai, tumor_cram, tumor_crai, intervals]
    }.set{cram_pair_intervals}


    no_intervals = false
    if (intervals == []) no_intervals = true

    manta_vcf            = Channel.empty()
    strelka_vcf          = Channel.empty()

    if (tools.contains('manta')) {
        MANTA(
            cram_pair,
            fasta,
            fasta_fai,
            target_bed_gz_tbi)

        manta_candidate_small_indels_vcf = MANTA.out.candidate_small_indels_vcf
        manta_candidate_sv_vcf           = MANTA.out.candidate_sv_vcf
        manta_diploid_sv_vcf             = MANTA.out.diploid_sv_vcf
        manta_somatic_sv_vcf             = MANTA.out.somatic_sv_vcf
        manta_csi_for_strelka_bp         = MANTA.out.manta_csi_for_strelka_bp

        manta_vcf = manta_candidate_small_indels_vcf.mix(manta_candidate_sv_vcf,manta_diploid_sv_vcf,manta_somatic_sv_vcf)

        if (tools.contains('strelka')) {
            STRELKA_BP(
                manta_csi_for_strelka_bp,
                fasta,
                fasta_fai,
                target_bed_gz_tbi)

            strelka_indels_vcf = STRELKA_BP.out.indels_vcf
            strelka_snvs_vcf   = STRELKA_BP.out.snvs_vcf

            strelka_vcf = strelka_vcf.mix(strelka_indels_vcf,strelka_snvs_vcf)
        }
    }

    if (tools.contains('msisensorpro')) {
        MSISENSORPRO_MSI_SOMATIC(
            cram_pair,
            msisensorpro_scan)
    }

    if (tools.contains('mutect2')){
        panel_of_normals.dump()
        MUTECT2(
            cram_pair_intervals,
            panel_of_normals,
            panel_of_normals_tbi,
            dict,
            fasta,
            fasta_fai,
            no_intervals,
            germline_resource,
            germline_resource_tbi
        )
    }

    if (tools.contains('strelka')) {
        STRELKA(
            cram_pair,
            fasta,
            fasta_fai,
            target_bed_gz_tbi)

        strelka_indels_vcf = STRELKA.out.indels_vcf
        strelka_snvs_vcf   = STRELKA.out.snvs_vcf

        strelka_vcf = strelka_vcf.mix(strelka_indels_vcf,strelka_snvs_vcf)
    }
    //tiddit should also work here

    emit:
        manta_vcf            = manta_vcf
        strelka_vcf          = strelka_vcf
}
