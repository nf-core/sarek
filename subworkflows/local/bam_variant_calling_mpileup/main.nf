include { CAT_CAT as CAT_MPILEUP           } from '../../../modules/nf-core/cat/cat/main'
include { BCFTOOLS_MPILEUP                 } from '../../../modules/nf-core/bcftools/mpileup/main'
include { GATK4_MERGEVCFS as MERGE_MPILEUP } from '../../../modules/nf-core/gatk4/mergevcfs/main'


workflow BAM_VARIANT_CALLING_MPILEUP {
    take:
        cram                    // channel: [mandatory] [meta, cram, interval]
        dict                    // channel: [mandatory]
        fasta                   // channel: [mandatory]

    main:

    ch_versions = Channel.empty()
    mpileup = Channel.empty()
    mpileup_vcf = Channel.empty()

    BCFTOOLS_MPILEUP(cram, fasta, true)
    BCFTOOLS_MPILEUP.out.mpileup.branch{
            intervals:    it[0].num_intervals > 1
            no_intervals: it[0].num_intervals <= 1
        }.set{mpileup_branch} 
    BCFTOOLS_MPILEUP.out.vcf.branch{
            intervals:    it[0].num_intervals > 1
            no_intervals: it[0].num_intervals <= 1
        }.set{mpileup_vcf_branch} 

    //Merge mpileup only when intervals and natural order sort them
    CAT_MPILEUP(mpileup_branch.intervals
        .map{ meta, pileup ->
            new_meta = meta.tumor_id ? [
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
                                            status:         meta.status,
                                            sex:            meta.sex,
                                        ]
            [groupKey(new_meta, meta.num_intervals), pileup]
            }
        .groupTuple(sort:true))

        MERGE_MPILEUP(
            mpileup_vcf_branch.intervals
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
            dict)
        mpileup_vcf = Channel.empty().mix(
            MERGE_MPILEUP.out.vcf,
            mpileup_vcf_branch.no_intervals)

    ch_versions = ch_versions.mix(BCFTOOLS_MPILEUP.out.versions)
    ch_versions = ch_versions.mix(CAT_MPILEUP.out.versions)
    ch_versions = ch_versions.mix(MERGE_MPILEUP.out.versions)

    emit:
    versions = ch_versions
    mpileup = Channel.empty().mix(CAT_MPILEUP.out.file_out, mpileup_branch.no_intervals)
    mpileup_vcf = mpileup_vcf
}
