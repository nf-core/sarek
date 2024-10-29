//
// PureCN variant calling
//
// For all modules here:
// A when clause condition is defined in the conf/modules.config to determine if the module should be run

include { GATK4_COLLECTREADCOUNTS } from '../../../modules/nf-core/gatk4/collectreadcounts/main'
include { GATK4_DENOISEREADCOUNTS } from '../../../modules/nf-core/gatk4/denoisereadcounts/main'
include { PURECN_RUN } from '../../../modules/nf-core/purecn/run/main'
include { PURECN_INTERVALFILE } from '../../../../nf-core-modules/modules/nf-core/purecn/intervalfile/main.nf'

workflow BAM_VARIANT_CALLING_SOMATIC_PURECN {

    take:
    tumor_file               // channel: [mandatory] [meta, tumor_cram, tumor_crai]
    fasta                    // channel: [mandatory] fasta needed for cram
    fasta_fai                // channel: [mandatory] fasta index needed by GATK
    dict                     // channel: [mandatory] fasta dictionary needed by GATK
    normaldb                 // channel: [mandatory] panel of normals built by PureCN
    gatk_pon                 // channel: [mandatory] panel of normals used by GATK for denoising
    intervals_bed            // channel: [mandatory] BED file processed by PureCN
    intervals_purecn         // channel: [mandatory] Interval file processed by PureCN

    main:

    ch_versions = Channel.empty()

    ch_fasta = [["id": "fasta"], fasta]
    ch_fai = [["id": "fasta_fai"], fasta_fai]
    ch_dict = [["id": "dict"], dict]
    ch_gatk_pon = [["id": "gatk_pon"], gatk_pon]

    ch_gatk_read_input = tumor_file.map {
        meta, tumor_cram, tumor_crai ->
        [meta, tumor_cram, tumor_crai, intervals_bed]
    }

    GATK4_COLLECTREADCOUNTS(ch_gatk_read_input, ch_fasta, ch_fai, ch_dict)
    ch_versions = ch_versions.mix(GATK4_COLLECTREADCOUNTS.out.versions)

    GATK4_DENOISEREADCOUNTS(GATK4_COLLECTREADCOUNTS.out.tsv, ch_gatk_pon)
    ch_versions = ch_versions.mix(GATK4_DENOISEREADCOUNTS.out.versions)

    ch_purecn_in = GATK4_DENOISEREADCOUNTS.out.denoised.map {
        meta, denoised -> [meta, denoised, intervals_purecn]
    }

    PURECN_RUN(ch_purecn_in, normaldb)
    ch_versions = ch_versions.mix(PURECN_RUN.out.versions)

    emit:
    versions = ch_versions
}
