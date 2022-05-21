include { TABIX_BGZIP as BGZIP_VC_STRELKA_INDELS } from '../../../../../modules/nf-core/modules/tabix/bgzip/main'
include { TABIX_BGZIP as BGZIP_VC_STRELKA_SNVS   } from '../../../../../modules/nf-core/modules/tabix/bgzip/main'
include { CONCAT_VCF as CONCAT_STRELKA_INDELS    } from '../../../../../modules/local/concat_vcf/main'
include { CONCAT_VCF as CONCAT_STRELKA_SNVS      } from '../../../../../modules/local/concat_vcf/main'
include { STRELKA_SOMATIC                        } from '../../../../../modules/nf-core/modules/strelka/somatic/main'

workflow RUN_STRELKA_SOMATIC {
    take:
    cram                     // channel: [mandatory] [meta, normal_cram, normal_crai, tumor_cram, tumor_crai, manta_vcf, manta_tbi, interval.bed.gz, interval.bed.gz.tbi] manta* are optional
    fasta                    // channel: [mandatory]
    fasta_fai                // channel: [mandatory]
    intervals_bed_gz         // channel: [optional]  Contains a bed.gz file of all intervals combined provided with the cram input(s). Mandatory if interval files are used.

    main:

    ch_versions = Channel.empty()

    STRELKA_SOMATIC(cram, fasta, fasta_fai )

    // Figure out if using intervals or no_intervals
    STRELKA_SOMATIC.out.vcf_snvs.branch{
            intervals:    it[0].num_intervals > 1
            no_intervals: it[0].num_intervals <= 1
        }.set{strelka_vcf_snvs}

    STRELKA_SOMATIC.out.vcf_indels.branch{
            intervals:    it[0].num_intervals > 1
            no_intervals: it[0].num_intervals <= 1
        }.set{strelka_vcf_indels}

    // Only when using intervals
    BGZIP_VC_STRELKA_SNVS(strelka_vcf_snvs.intervals)

    CONCAT_STRELKA_SNVS(BGZIP_VC_STRELKA_SNVS.out.output.map{ meta, vcf ->
                new_meta = [patient:meta.patient, normal_id:meta.normal_id, tumor_id:meta.tumor_id, gender:meta.gender, id:meta.tumor_id + "_vs_" + meta.normal_id, num_intervals:meta.num_intervals]

                [groupKey(new_meta, new_meta.num_intervals), vcf]
            }.groupTuple(),
            fasta_fai,
            intervals_bed_gz)

    BGZIP_VC_STRELKA_INDELS(strelka_vcf_indels.intervals)

    CONCAT_STRELKA_INDELS(BGZIP_VC_STRELKA_INDELS.out.output.map{ meta, vcf ->
                new_meta = [patient:meta.patient, normal_id:meta.normal_id, tumor_id:meta.tumor_id, gender:meta.gender, id:meta.tumor_id + "_vs_" + meta.normal_id, num_intervals:meta.num_intervals]

                [groupKey(new_meta, new_meta.num_intervals), vcf]
            }.groupTuple(),
            fasta_fai,
            intervals_bed_gz)

    // Mix output channels for "no intervals" and "with intervals" results
    strelka_vcf = Channel.empty().mix(
                    CONCAT_STRELKA_SNVS.out.vcf,
                    CONCAT_STRELKA_INDELS.out.vcf,
                    strelka_vcf_snvs.no_intervals,
                    strelka_vcf_indels.no_intervals)
                .map{ meta, vcf ->
                    [[patient:meta.patient, normal_id:meta.normal_id, tumor_id:meta.tumor_id, gender:meta.gender, id:meta.tumor_id + "_vs_" + meta.normal_id, num_intervals:meta.num_intervals, variantcaller:"Strelka"]
                    , vcf]
                }

    ch_versions = ch_versions.mix(BGZIP_VC_STRELKA_SNVS.out.versions)
    ch_versions = ch_versions.mix(BGZIP_VC_STRELKA_INDELS.out.versions)
    ch_versions = ch_versions.mix(CONCAT_STRELKA_SNVS.out.versions)
    ch_versions = ch_versions.mix(CONCAT_STRELKA_INDELS.out.versions)
    ch_versions = ch_versions.mix(STRELKA_SOMATIC.out.versions)

    emit:
    strelka_vcf
    versions = ch_versions
}
