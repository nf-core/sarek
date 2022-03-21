include { BGZIP as BGZIP_VC_MANTA_SMALL_INDELS    } from '../../../modules/local/bgzip'
include { BGZIP as BGZIP_VC_MANTA_SV              } from '../../../modules/local/bgzip'
include { BGZIP as BGZIP_VC_MANTA_TUMOR           } from '../../../modules/local/bgzip'
include { CONCAT_VCF as CONCAT_MANTA_SMALL_INDELS } from '../../../modules/local/concat_vcf/main'
include { CONCAT_VCF as CONCAT_MANTA_SV           } from '../../../modules/local/concat_vcf/main'
include { CONCAT_VCF as CONCAT_MANTA_TUMOR        } from '../../../modules/local/concat_vcf/main'
include { MANTA_TUMORONLY                         } from '../../../modules/local/manta/tumoronly/main'

workflow RUN_MANTA_TUMORONLY {
    take:
    cram_recalibrated_intervals_gz_tbi
    fasta
    fasta_fai
    num_intervals
    intervals_bed_combine_gz

    main:

    ch_versions = Channel.empty()
    MANTA_TUMORONLY(
            cram_recalibrated_intervals_gz_tbi,
            fasta,
            fasta_fai
        )

        ch_versions = ch_versions.mix(MANTA_TUMORONLY.out.versions)

        if(no_intervals){
            manta_candidate_small_indels_vcf = MANTA_TUMORONLY.out.candidate_small_indels_vcf
            manta_candidate_sv_vcf           = MANTA_TUMORONLY.out.candidate_sv_vcf
            manta_tumor_sv_vcf               = MANTA_TUMORONLY.out.tumor_sv_vcf
        }else{

            BGZIP_VC_MANTA_SV(MANTA_TUMORONLY.out.candidate_small_indels_vcf)
            BGZIP_VC_MANTA_SMALL_INDELS(MANTA_TUMORONLY.out.candidate_sv_vcf)
            BGZIP_VC_MANTA_TUMOR(MANTA_TUMORONLY.out.tumor_sv_vcf)

            BGZIP_VC_MANTA_SV.out.vcf.map{ meta, vcf ->
                new_meta = meta.clone()
                new_meta.id = new_meta.sample
                [new_meta, vcf]
            }.groupTuple(size: num_intervals)
            .set{manta_sv_vcf_to_concat}

            BGZIP_VC_MANTA_SMALL_INDELS.out.vcf.map{ meta, vcf ->
                new_meta = meta.clone()
                new_meta.id = new_meta.sample
                [new_meta, vcf]
            }.groupTuple(size: num_intervals)
            .set{manta_small_indels_vcf_to_concat}

            BGZIP_VC_MANTA_TUMOR.out.vcf.map{ meta, vcf ->
                new_meta = meta.clone()
                new_meta.id = new_meta.sample
                [new_meta, vcf]
            }.groupTuple(size: num_intervals)
            .set{manta_tumor_sv_vcf_to_concat}

            CONCAT_MANTA_SV(manta_sv_vcf_to_concat, fasta_fai, intervals_bed_combine_gz)
            CONCAT_MANTA_SMALL_INDELS(manta_small_indels_vcf_to_concat,fasta_fai, intervals_bed_combine_gz)
            CONCAT_MANTA_TUMOR(manta_tumor_sv_vcf_to_concat, fasta_fai, intervals_bed_combine_gz)

            manta_candidate_small_indels_vcf = CONCAT_MANTA_SV.out.vcf
            manta_candidate_sv_vcf           = CONCAT_MANTA_SMALL_INDELS.out.vcf
            manta_tumor_sv_vcf               = CONCAT_MANTA_TUMOR.out.vcf

            ch_versions = ch_versions.mix(BGZIP_VC_MANTA_SV.out.versions)
            ch_versions = ch_versions.mix(BGZIP_VC_MANTA_SMALL_INDELS.out.versions)
            ch_versions = ch_versions.mix(BGZIP_VC_MANTA_TUMOR.out.versions)

            ch_versions = ch_versions.mix(CONCAT_MANTA_SV.out.versions)
            ch_versions = ch_versions.mix(CONCAT_MANTA_SMALL_INDELS.out.versions)
            ch_versions = ch_versions.mix(CONCAT_MANTA_TUMOR.out.versions)
        }

        manta_vcf = manta_vcf.mix(manta_candidate_small_indels_vcf, manta_candidate_sv_vcf, manta_tumor_sv_vcf)

    emit:
    versions = ch_versions
    manta_vcf
}
