include { TABIX_BGZIP as BGZIP_VC_MANTA_DIPLOID      } from '../../../../../modules/nf-core/modules/tabix/bgzip/main'
include { TABIX_BGZIP as BGZIP_VC_MANTA_SMALL_INDELS } from '../../../../../modules/nf-core/modules/tabix/bgzip/main'
include { TABIX_BGZIP as BGZIP_VC_MANTA_SOMATIC      } from '../../../../../modules/nf-core/modules/tabix/bgzip/main'
include { TABIX_BGZIP as BGZIP_VC_MANTA_SV           } from '../../../../../modules/nf-core/modules/tabix/bgzip/main'
include { CONCAT_VCF as CONCAT_MANTA_DIPLOID         } from '../../../../../modules/local/concat_vcf/main'
include { CONCAT_VCF as CONCAT_MANTA_SMALL_INDELS    } from '../../../../../modules/local/concat_vcf/main'
include { CONCAT_VCF as CONCAT_MANTA_SOMATIC         } from '../../../../../modules/local/concat_vcf/main'
include { CONCAT_VCF as CONCAT_MANTA_SV              } from '../../../../../modules/local/concat_vcf/main'
include { MANTA_SOMATIC                              } from '../../../../../modules/nf-core/modules/manta/somatic/main'

workflow RUN_MANTA_SOMATIC {
    take:
    cram                     // channel: [mandatory] [meta, normal_cram, normal_crai, tumor_cram, tumor_crai, interval.bed.gz, interval.bed.gz.tbi]
    fasta                    // channel: [mandatory]
    fasta_fai                // channel: [mandatory]
    intervals_bed_gz         // channel: [optional]  Contains a bed.gz file of all intervals combined provided with the cram input(s). Mandatory if interval files are used.
    num_intervals            //     val: [optional]  Number of used intervals, mandatory when intervals are provided.

    main:

    ch_versions = Channel.empty()

    MANTA_SOMATIC(cram, fasta, fasta_fai)

    // Figure out if using intervals or no_intervals
    MANTA_SOMATIC.out.candidate_small_indels_vcf.branch{
            intervals:    num_intervals > 1
            no_intervals: num_intervals == 1
        }.set{manta_candidate_small_indels_vcf}

    MANTA_SOMATIC.out.candidate_small_indels_vcf_tbi.branch{
            intervals:    num_intervals > 1
            no_intervals: num_intervals == 1
        }.set{manta_candidate_small_indels_vcf_tbi}

    MANTA_SOMATIC.out.candidate_sv_vcf.branch{
            intervals:    num_intervals > 1
            no_intervals: num_intervals == 1
        }.set{manta_candidate_sv_vcf}

    MANTA_SOMATIC.out.diploid_sv_vcf.branch{
            intervals:    num_intervals > 1
            no_intervals: num_intervals == 1
        }.set{manta_diploid_sv_vcf}

    MANTA_SOMATIC.out.somatic_sv_vcf.branch{
            intervals:    num_intervals > 1
            no_intervals: num_intervals == 1
        }.set{manta_somatic_sv_vcf}

    //Only when using intervals
    BGZIP_VC_MANTA_SV(manta_candidate_small_indels_vcf.intervals)

    CONCAT_MANTA_SV(
        BGZIP_VC_MANTA_SV.out.output.map{ meta, vcf ->
                new_meta = meta.clone()
                new_meta.id = new_meta.tumor_id + "_vs_" + new_meta.normal_id
                [new_meta, vcf]
            }.groupTuple(size: num_intervals),
        fasta_fai,
        intervals_bed_gz)

    BGZIP_VC_MANTA_SMALL_INDELS(manta_candidate_sv_vcf.intervals)

    CONCAT_MANTA_SMALL_INDELS(
        BGZIP_VC_MANTA_SMALL_INDELS.out.output.map{ meta, vcf ->
                new_meta = meta.clone()
                new_meta.id = new_meta.tumor_id + "_vs_" + new_meta.normal_id
                [new_meta, vcf]
            }.groupTuple(size: num_intervals),
        fasta_fai,
        intervals_bed_gz)

    BGZIP_VC_MANTA_DIPLOID(manta_diploid_sv_vcf.intervals)

    CONCAT_MANTA_DIPLOID(
        BGZIP_VC_MANTA_DIPLOID.out.output.map{ meta, vcf ->
                new_meta = meta.clone()
                new_meta.id = new_meta.tumor_id + "_vs_" + new_meta.normal_id
                [new_meta, vcf]
            }.groupTuple(size: num_intervals),
        fasta_fai,
        intervals_bed_gz)

    BGZIP_VC_MANTA_SOMATIC(manta_somatic_sv_vcf.intervals)

    CONCAT_MANTA_SOMATIC(
        BGZIP_VC_MANTA_SOMATIC.out.output.map{ meta, vcf ->
                new_meta = meta.clone()
                new_meta.id = new_meta.tumor_id + "_vs_" + new_meta.normal_id
                [new_meta, vcf]
            }.groupTuple(size: num_intervals),
        fasta_fai,
        intervals_bed_gz)

    // Mix output channels for "no intervals" and "with intervals" results
    manta_vcf = Channel.empty().mix(
        CONCAT_MANTA_SV.out.vcf,
        CONCAT_MANTA_SMALL_INDELS.out.vcf,
        CONCAT_MANTA_DIPLOID.out.vcf,
        CONCAT_MANTA_SOMATIC.out.vcf,
        manta_candidate_sv_vcf.no_intervals,
        manta_candidate_small_indels_vcf.no_intervals,
        manta_diploid_sv_vcf.no_intervals,
        manta_somatic_sv_vcf.no_intervals
    ).map{ meta, vcf ->
        meta.variantcaller = "Manta"
        [meta, vcf]
    }

    manta_candidate_small_indels_vcf = Channel.empty().mix(
        CONCAT_MANTA_SMALL_INDELS.out.vcf,
        manta_candidate_small_indels_vcf.no_intervals
    ).map{ meta, vcf ->
        meta.variantcaller = "Manta"
        [meta, vcf]
    }

    manta_candidate_small_indels_vcf_tbi = Channel.empty().mix(
        CONCAT_MANTA_SMALL_INDELS.out.tbi,
        manta_candidate_small_indels_vcf_tbi.no_intervals
    ).map{ meta, vcf ->
        meta.variantcaller = "Manta"
        [meta, vcf]
    }

    ch_versions = ch_versions.mix(BGZIP_VC_MANTA_SV.out.versions)
    ch_versions = ch_versions.mix(BGZIP_VC_MANTA_SMALL_INDELS.out.versions)
    ch_versions = ch_versions.mix(BGZIP_VC_MANTA_DIPLOID.out.versions)
    ch_versions = ch_versions.mix(BGZIP_VC_MANTA_SOMATIC.out.versions)
    ch_versions = ch_versions.mix(CONCAT_MANTA_SV.out.versions)
    ch_versions = ch_versions.mix(CONCAT_MANTA_SMALL_INDELS.out.versions)
    ch_versions = ch_versions.mix(CONCAT_MANTA_DIPLOID.out.versions)
    ch_versions = ch_versions.mix(CONCAT_MANTA_SOMATIC.out.versions)
    ch_versions = ch_versions.mix(MANTA_SOMATIC.out.versions)

    emit:
    manta_vcf
    manta_candidate_small_indels_vcf
    manta_candidate_small_indels_vcf_tbi
    versions = ch_versions

}
