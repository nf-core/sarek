include { BGZIP as BGZIP_VC_MANTA_DIPLOID           } from '../../../modules/local/bgzip'
include { BGZIP as BGZIP_VC_MANTA_SMALL_INDELS      } from '../../../modules/local/bgzip'
include { BGZIP as BGZIP_VC_MANTA_SOMATIC           } from '../../../modules/local/bgzip'
include { BGZIP as BGZIP_VC_MANTA_SV                } from '../../../modules/local/bgzip'
include { CONCAT_VCF as CONCAT_MANTA_DIPLOID        } from '../../../modules/local/concat_vcf/main'
include { CONCAT_VCF as CONCAT_MANTA_SMALL_INDELS   } from '../../../modules/local/concat_vcf/main'
include { CONCAT_VCF as CONCAT_MANTA_SOMATIC        } from '../../../modules/local/concat_vcf/main'
include { CONCAT_VCF as CONCAT_MANTA_SV             } from '../../../modules/local/concat_vcf/main'
include { MANTA_SOMATIC                             } from '../../../modules/local/manta/somatic/main'

workflow RUN_MANTA_SOMATIC {
    take:
    cram_pair_intervals_gz_tbi
    fasta
    fasta_fai
    num_intervals
    intervals_bed_combine_gz

    main:

    ch_versions = Channel.empty()
        MANTA_SOMATIC(
            cram_pair_intervals_gz_tbi,
            fasta,
            fasta_fai)

        ch_versions = ch_versions.mix(MANTA_SOMATIC.out.versions)

        if (no_intervals) {
            manta_candidate_small_indels_vcf = MANTA_SOMATIC.out.candidate_small_indels_vcf
            manta_candidate_sv_vcf           = MANTA_SOMATIC.out.candidate_sv_vcf
            manta_diploid_sv_vcf             = MANTA_SOMATIC.out.diploid_sv_vcf
            manta_somatic_sv_vcf             = MANTA_SOMATIC.out.somatic_sv_vcf
        } else {
            BGZIP_VC_MANTA_SV(MANTA_SOMATIC.out.candidate_small_indels_vcf)
            BGZIP_VC_MANTA_SMALL_INDELS(MANTA_SOMATIC.out.candidate_sv_vcf)
            BGZIP_VC_MANTA_DIPLOID(MANTA_SOMATIC.out.diploid_sv_vcf)
            BGZIP_VC_MANTA_SOMATIC(MANTA_SOMATIC.out.somatic_sv_vcf)

            manta_sv_vcf_to_concat = BGZIP_VC_MANTA_SV.out.vcf.map{ meta, vcf ->
                new_meta = meta.clone()
                new_meta.id = new_meta.tumor_id + "_vs_" + new_meta.normal_id
                [new_meta, vcf]
            }.groupTuple(size: num_intervals)

            manta_small_indels_vcf_to_concat = BGZIP_VC_MANTA_SMALL_INDELS.out.vcf.map{ meta, vcf ->
                new_meta = meta.clone()
                new_meta.id = new_meta.tumor_id + "_vs_" + new_meta.normal_id
                [new_meta, vcf]
            }.groupTuple(size: num_intervals)

            manta_diploid_vcf_to_concat = BGZIP_VC_MANTA_DIPLOID.out.vcf.map{ meta, vcf ->
                new_meta = meta.clone()
                new_meta.id = new_meta.tumor_id + "_vs_" + new_meta.normal_id
                [new_meta, vcf]
            }.groupTuple(size: num_intervals)

            manta_somatic_sv_vcf_to_concat = BGZIP_VC_MANTA_SOMATIC.out.vcf.map{ meta, vcf ->
                new_meta = meta.clone()
                new_meta.id = new_meta.tumor_id + "_vs_" + new_meta.normal_id
                [new_meta, vcf]
            }.groupTuple(size: num_intervals)

            CONCAT_MANTA_SV(manta_sv_vcf_to_concat, fasta_fai, intervals_bed_combine_gz)
            CONCAT_MANTA_SMALL_INDELS(manta_small_indels_vcf_to_concat,fasta_fai, intervals_bed_combine_gz)
            CONCAT_MANTA_DIPLOID(manta_diploid_vcf_to_concat, fasta_fai, intervals_bed_combine_gz)
            CONCAT_MANTA_SOMATIC(manta_somatic_sv_vcf_to_concat, fasta_fai, intervals_bed_combine_gz)

            manta_candidate_small_indels_vcf = CONCAT_MANTA_SV.out.vcf
            manta_candidate_sv_vcf           = CONCAT_MANTA_SMALL_INDELS.out.vcf
            manta_diploid_sv_vcf             = CONCAT_MANTA_DIPLOID.out.vcf
            manta_somatic_sv_vcf             = CONCAT_MANTA_SOMATIC.out.vcf

            ch_versions = ch_versions.mix(BGZIP_VC_MANTA_SV.out.versions)
            ch_versions = ch_versions.mix(BGZIP_VC_MANTA_SMALL_INDELS.out.versions)
            ch_versions = ch_versions.mix(BGZIP_VC_MANTA_DIPLOID.out.versions)
            ch_versions = ch_versions.mix(BGZIP_VC_MANTA_SOMATIC.out.versions)

            ch_versions = ch_versions.mix(CONCAT_MANTA_SV.out.versions)
            ch_versions = ch_versions.mix(CONCAT_MANTA_SMALL_INDELS.out.versions)
            ch_versions = ch_versions.mix(CONCAT_MANTA_DIPLOID.out.versions)
            ch_versions = ch_versions.mix(CONCAT_MANTA_SOMATIC.out.versions)

        }

        manta_vcf = manta_vcf.mix(manta_candidate_small_indels_vcf,manta_candidate_sv_vcf,manta_diploid_sv_vcf,manta_somatic_sv_vcf)

    emit:
    versions = ch_versions

}
