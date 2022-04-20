include { BGZIP as BGZIP_VC_MANTA_SMALL_INDELS    } from '../../../../../modules/nf-core/modules/tabix/bgzip/main'
include { BGZIP as BGZIP_VC_MANTA_SV              } from '../../../../../modules/nf-core/modules/tabix/bgzip/main'
include { BGZIP as BGZIP_VC_MANTA_TUMOR           } from '../../../../../modules/nf-core/modules/tabix/bgzip/main'
include { CONCAT_VCF as CONCAT_MANTA_SMALL_INDELS } from '../../../../../modules/local/concat_vcf/main'
include { CONCAT_VCF as CONCAT_MANTA_SV           } from '../../../../../modules/local/concat_vcf/main'
include { CONCAT_VCF as CONCAT_MANTA_TUMOR        } from '../../../../../modules/local/concat_vcf/main'
include { MANTA_TUMORONLY                         } from '../../../../../modules/nf-core/modules/manta/tumoronly/main'

// TODO: Research if splitting by intervals is ok, we pretend for now it is fine.
// Seems to be the consensus on upstream modules implementation too
workflow RUN_MANTA_TUMORONLY {
    take:
    cram                     // channel: [mandatory] [meta, cram, crai, interval.bed.gz, interval.bed.gz.tbi]
    fasta                    // channel: [mandatory]
    fasta_fai                // channel: [mandatory]
    intervals_bed_gz         // channel: [optional]  Contains a bed.gz file of all intervals combined provided with the cram input(s). Mandatory if interval files are used.
    num_intervals            //     val: [optional]  Number of used intervals, mandatory when intervals are provided.

    main:

    ch_versions = Channel.empty()

    MANTA_TUMORONLY(cram, fasta, fasta_fai)

    // Figure out if using intervals or no_intervals
    MANTA_TUMORONLY.out.candidate_small_indels_vcf.branch{
            intervals:    num_intervals > 1
            no_intervals: num_intervals == 1
        }.set{manta_small_indels_vcf}

    MANTA_TUMORONLY.out.candidate_sv_vcf.branch{
            intervals:    num_intervals > 1
            no_intervals: num_intervals == 1
        }.set{manta_candidate_sv_vcf}

    MANTA_TUMORONLY.out.tumor_sv_vcf.branch{
            intervals:    num_intervals > 1
            no_intervals: num_intervals == 1
        }.set{manta_tumor_sv_vcf}

    //Only when using intervals
    BGZIP_VC_MANTA_SMALL_INDELS(manta_small_indels_vcf.intervals)

    CONCAT_MANTA_SMALL_INDELS(
        BGZIP_VC_MANTA_SMALL_INDELS.out.output.map{ meta, vcf ->
                new_meta = meta.clone()
                new_meta.id = new_meta.sample
                [new_meta, vcf]
            }.groupTuple(size: num_intervals),
        fasta_fai,
        intervals_bed_gz)

    BGZIP_VC_MANTA_SV(manta_candidate_sv_vcf.intervals)

    CONCAT_MANTA_SV(
        BGZIP_VC_MANTA_SV.out.output.map{ meta, vcf ->
                new_meta = meta.clone()
                new_meta.id = new_meta.sample
                [new_meta, vcf]
            }.groupTuple(size: num_intervals),
        fasta_fai,
        intervals_bed_gz)

    BGZIP_VC_MANTA_TUMOR(manta_tumor_sv_vcf.intervals)

    CONCAT_MANTA_TUMOR(
        BGZIP_VC_MANTA_TUMOR.out.output.map{ meta, vcf ->
                new_meta = meta.clone()
                new_meta.id = new_meta.sample
                [new_meta, vcf]
            }.groupTuple(size: num_intervals),
        fasta_fai,
        intervals_bed_gz)

    // Mix output channels for "no intervals" and "with intervals" results
    manta_vcf = Channel.empty().mix(
        CONCAT_MANTA_SMALL_INDELS.out.vcf,
        CONCAT_MANTA_SV.out.vcf,
        CONCAT_MANTA_TUMOR.out.vcf,
        manta_small_indels_vcf.no_intervals,
        manta_candidate_sv_vcf.no_intervals,
        manta_tumor_sv_vcf.no_intervals
    )

    ch_versions = ch_versions.mix(BGZIP_VC_MANTA_SV.out.versions)
    ch_versions = ch_versions.mix(BGZIP_VC_MANTA_SMALL_INDELS.out.versions)
    ch_versions = ch_versions.mix(BGZIP_VC_MANTA_TUMOR.out.versions)
    ch_versions = ch_versions.mix(CONCAT_MANTA_SV.out.versions)
    ch_versions = ch_versions.mix(CONCAT_MANTA_SMALL_INDELS.out.versions)
    ch_versions = ch_versions.mix(CONCAT_MANTA_TUMOR.out.versions)
    ch_versions = ch_versions.mix(MANTA_TUMORONLY.out.versions)

    emit:
    manta_vcf
    versions = ch_versions
}
