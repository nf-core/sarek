include { CAT_CAT as CAT_MPILEUP         } from '../../../modules/nf-core/cat/cat/main'
include { BCFTOOLS_MPILEUP               } from '../../../modules/nf-core/bcftools/mpileup/main'
include { SAMTOOLS_MPILEUP               } from '../../../modules/nf-core/samtools/mpileup/main'
include { GATK4_MERGEVCFS                } from '../../../modules/nf-core/gatk4/mergevcfs/main'


workflow BAM_VARIANT_CALLING_MPILEUP {
    take:
        cram                    // channel: [mandatory] [meta, cram, interval]
        fasta                   // channel: [mandatory]
        dict

    main:

    ch_versions = Channel.empty()

    keep_bcftools_mpileup = false
    BCFTOOLS_MPILEUP(cram, fasta, keep_bcftools_mpileup)
    SAMTOOLS_MPILEUP(cram, fasta)
    vcfs = BCFTOOLS_MPILEUP.out.vcf.branch{
            intervals:    it[0].num_intervals > 1
            no_intervals: it[0].num_intervals <= 1
        }
    mpileup = SAMTOOLS_MPILEUP.out.mpileup.branch{
            intervals:    it[0].num_intervals > 1
            no_intervals: it[0].num_intervals <= 1
        }

    //Merge mpileup only when intervals and natural order sort them
    CAT_MPILEUP(mpileup.intervals
        .map{ meta, pileup ->
            def new_meta = meta.tumor_id ? [
                                            id:             meta.tumor_id + "_vs_" + meta.normal_id,
                                            normal_id:      meta.normal_id,
                                            num_intervals:  meta.num_intervals,
                                            patient:        meta.patient,
                                            sex:            meta.sex,
                                            tumor_id:       meta.tumor_id,
                                        ] // not annotated, so no variantcaller necessary
                                        : [
                                            id:             meta.sample,
                                            num_intervals:  meta.num_intervals,
                                            patient:        meta.patient,
                                            sample:         meta.sample,
                                            sex:            meta.sex,
                                            status:         meta.status
                                        ]
            [groupKey(new_meta, meta.num_intervals), pileup]
            }
        .groupTuple(sort:true))

    GATK4_MERGEVCFS(vcfs.intervals
        .map{ meta, vcf ->
            new_meta = [
                        id:             meta.sample,
                        num_intervals:  meta.num_intervals,
                        patient:        meta.patient,
                        sample:         meta.sample,
                        sex:            meta.sex,
                        status:         meta.status
                    ]

            [groupKey(new_meta, new_meta.num_intervals), vcf]
        }.groupTuple(),
    dict.map{ it -> [[id:it[0].baseName], it]})

    ch_versions = ch_versions.mix(SAMTOOLS_MPILEUP.out.versions)
    ch_versions = ch_versions.mix(BCFTOOLS_MPILEUP.out.versions)
    ch_versions = ch_versions.mix(CAT_MPILEUP.out.versions)
    ch_versions = ch_versions.mix(GATK4_MERGEVCFS.out.versions)

    emit:
    versions = ch_versions
    mpileup = Channel.empty().mix(CAT_MPILEUP.out.file_out, mpileup.no_intervals)
    vcf = Channel.empty().mix(GATK4_MERGEVCFS.out.vcf, vcfs.no_intervals)
}
