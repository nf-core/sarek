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
    versions = Channel.empty()

    STRELKA_SOMATIC(cram, fasta, fasta_fai )

    // Figure out if using intervals or no_intervals
    vcf_snvs = STRELKA_SOMATIC.out.vcf_snvs.branch{
        intervals:    it[0].num_intervals > 1
        no_intervals: it[0].num_intervals <= 1
    }

    vcf_indels = STRELKA_SOMATIC.out.vcf_indels.branch{
        intervals:    it[0].num_intervals > 1
        no_intervals: it[0].num_intervals <= 1
    }

    // Only when using intervals
    MERGE_STRELKA_SNVS(vcf_snvs.intervals.map{ meta, vcf ->

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
            dict.map{ it -> [ [ id:'dict' ], it ] })

    MERGE_STRELKA_INDELS(vcf_indels.intervals.map{ meta, vcf ->

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
            dict.map{ it -> [ [ id:'dict' ], it ] })

    // Mix output channels for "no intervals" and "with intervals" results
    vcf = Channel.empty().mix(
                    MERGE_STRELKA_SNVS.out.vcf,
                    vcf_snvs.no_intervals,
                    MERGE_STRELKA_INDELS.out.vcf,
                    vcf_indels.no_intervals
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

    versions = versions.mix(MERGE_STRELKA_SNVS.out.versions)
    versions = versions.mix(MERGE_STRELKA_INDELS.out.versions)
    versions = versions.mix(STRELKA_SOMATIC.out.versions)

    emit:
    vcf

    versions
}
