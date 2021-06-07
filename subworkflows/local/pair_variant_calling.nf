/*
========================================================================================
    SOMATIC VARIANT CALLING
========================================================================================
*/

params.manta_options                  = [:]
params.msisensorpro_msi_options       = [:]
params.strelka_options                = [:]
params.strelka_bp_options             = [:]

include { MANTA_SOMATIC as MANTA }                       from '../../modules/nf-core/software/manta/somatic/main'     addParams(options: params.manta_options)
include { MSISENSORPRO_MSI }                             from '../../modules/nf-core/software/msisensorpro/msi/main'  addParams(options: params.msisensorpro_msi_options)
include { STRELKA_SOMATIC as STRELKA }                   from '../../modules/nf-core/software/strelka/somatic/main'   addParams(options: params.strelka_options)
include { STRELKA_SOMATIC_BEST_PRACTICES as STRELKA_BP } from '../../modules/nf-core/software/strelka/somaticbp/main' addParams(options: params.strelka_bp_options)

workflow PAIR_VARIANT_CALLING {
    take:
        bam               // channel: [mandatory] bam
        dbsnp             // channel: [mandatory] dbsnp
        dbsnp_tbi         // channel: [mandatory] dbsnp_tbi
        dict              // channel: [mandatory] dict
        fai               // channel: [mandatory] fai
        fasta             // channel: [mandatory] fasta
        intervals         // channel: [mandatory] intervals
        msisensor_scan    // channel: [optional]  msisensor_scan
        target_bed        // channel: [optional]  target_bed
        target_bed_gz_tbi // channel: [optional]  target_bed_gz_tbi

    main:

    bam.map{ meta, bam, bai ->
        patient = meta.patient
        sample  = meta.sample
        gender  = meta.gender
        status  = meta.status
        [patient, sample, gender, status, bam, bai]
    }.branch{
        normal: it[3] == 0
        tumor:  it[3] == 1
    }.set{ bam_to_cross }

    bam_pair = bam_to_cross.normal.cross(bam_to_cross.tumor).map { normal, tumor ->
        def meta = [:]
        meta.patient = normal[0]
        meta.normal  = normal[1]
        meta.tumor   = tumor[1]
        meta.gender  = normal[2]
        meta.id      = "${meta.tumor}_vs_${meta.normal}".toString()

        [meta, normal[4], normal[5], tumor[4], tumor[5]]
    }

    manta_vcf            = Channel.empty()
    strelka_vcf          = Channel.empty()

    if ('manta' in params.tools.toLowerCase()) {
        MANTA(
            bam_pair,
            fasta,
            fai,
            target_bed_gz_tbi)

        manta_candidate_small_indels_vcf = MANTA.out.candidate_small_indels_vcf
        manta_candidate_sv_vcf           = MANTA.out.candidate_sv_vcf
        manta_diploid_sv_vcf             = MANTA.out.diploid_sv_vcf
        manta_somatic_sv_vcf             = MANTA.out.somatic_sv_vcf
        manta_csi_for_strelka_bp         = MANTA.out.manta_csi_for_strelka_bp

        manta_vcf = manta_candidate_small_indels_vcf.mix(manta_candidate_sv_vcf,manta_diploid_sv_vcf,manta_somatic_sv_vcf)

        if ('strelka' in params.tools.toLowerCase()) {
            STRELKA_BP(
                manta_csi_for_strelka_bp,
                fasta,
                fai,
                target_bed_gz_tbi)

            strelka_indels_vcf = STRELKA_BP.out.indels_vcf
            strelka_snvs_vcf   = STRELKA_BP.out.snvs_vcf

            strelka_vcf = strelka_vcf.mix(strelka_indels_vcf,strelka_snvs_vcf)
        }
    }

    if ('msisensor' in params.tools.toLowerCase()) {
        MSISENSORPRO_MSI(
            bam_pair,
            msisensor_scan)
    }

    if ('strelka' in params.tools.toLowerCase()) {
        STRELKA(
            bam_pair,
            fasta,
            fai,
            target_bed_gz_tbi)

        strelka_indels_vcf = STRELKA.out.indels_vcf
        strelka_snvs_vcf   = STRELKA.out.snvs_vcf

        strelka_vcf = strelka_vcf.mix(strelka_indels_vcf,strelka_snvs_vcf)
    }

    emit:
        manta_vcf            = manta_vcf
        strelka_vcf          = strelka_vcf
}
