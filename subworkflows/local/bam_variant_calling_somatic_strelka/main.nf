include { GATK4_MERGEVCFS as MERGE_STRELKA_INDELS } from '../../../modules/nf-core/gatk4/mergevcfs/main'
include { GATK4_MERGEVCFS as MERGE_STRELKA_SNVS   } from '../../../modules/nf-core/gatk4/mergevcfs/main'
include { STRELKA_SOMATIC                         } from '../../../modules/nf-core/strelka/somatic/main'

workflow BAM_VARIANT_CALLING_SOMATIC_STRELKA {
    take:
    cram                     // channel: [mandatory] [meta, normal_cram, normal_crai, tumor_cram, tumor_crai, manta_vcf, manta_tbi, interval.bed.gz, interval.bed.gz.tbi] manta* are optional
    dict                     // channel: [optional]
    fasta                    // channel: [mandatory]
    fasta_fai                // channel: [mandatory]

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
    MERGE_STRELKA_SNVS(strelka_vcf_snvs.intervals.map{ meta, vcf ->

                [groupKey([
                            id:             meta.tumor_id + "_vs_" + meta.normal_id,
                            normal_id:      meta.normal_id,
                            num_intervals:  meta.num_intervals,
                            patient:        meta.patient,
                            sex:            meta.sex,
                            tumor_id:       meta.tumor_id,
                            ],
                        meta.num_intervals),
                vcf]

            }.groupTuple(),
            dict)

    MERGE_STRELKA_INDELS(strelka_vcf_indels.intervals.map{ meta, vcf ->

                [groupKey([
                            id:             meta.tumor_id + "_vs_" + meta.normal_id,
                            normal_id:      meta.normal_id,
                            num_intervals:  meta.num_intervals,
                            patient:        meta.patient,
                            sex:            meta.sex,
                            tumor_id:       meta.tumor_id,
                            ],
                            meta.num_intervals),
                vcf]
            }.groupTuple(),
            dict)

    // Mix output channels for "no intervals" and "with intervals" results
    strelka_vcf = Channel.empty().mix(
                    MERGE_STRELKA_SNVS.out.vcf,
                    strelka_vcf_snvs.no_intervals,
                    MERGE_STRELKA_INDELS.out.vcf,
                    strelka_vcf_indels.no_intervals
                    )
                .map{ meta, vcf ->
                    [[
                        id:             meta.tumor_id + "_vs_" + meta.normal_id,
                        normal_id:      meta.normal_id,
                        num_intervals:  meta.num_intervals,
                        patient:        meta.patient,
                        sex:            meta.sex,
                        tumor_id:       meta.tumor_id,
                        variantcaller:  "strelka"
                        ],
                    vcf]
                }

    ch_versions = ch_versions.mix(MERGE_STRELKA_SNVS.out.versions)
    ch_versions = ch_versions.mix(MERGE_STRELKA_INDELS.out.versions)
    ch_versions = ch_versions.mix(STRELKA_SOMATIC.out.versions)

    emit:
    strelka_vcf
    versions = ch_versions
}
