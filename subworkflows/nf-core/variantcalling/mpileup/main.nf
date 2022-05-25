include { CAT_CAT as CAT_MPILEUP_NORMAL                   } from '../../../../modules/nf-core/modules/cat/cat/main.nf'
include { CAT_CAT as CAT_MPILEUP_TUMOR                    } from '../../../../modules/nf-core/modules/cat/cat/main.nf'
include { SAMTOOLS_MPILEUP as MPILEUP_NORMAL              } from '../../../../modules/nf-core/modules/samtools/mpileup/main'
include { SAMTOOLS_MPILEUP as MPILEUP_TUMOR               } from '../../../../modules/nf-core/modules/samtools/mpileup/main'

workflow MPILEUP {
    take:
    cram_normal              // channel: [mandatory] [meta, cram, crai, interval]
    cram_tumor               // channel: [mandatory]
    fasta                    // channel: [mandatory]

    main:

    ch_versions = Channel.empty()

    MPILEUP_NORMAL(cram_normal, fasta)

    MPILEUP_NORMAL.out.mpileup.branch{
            intervals:    it[0].num_intervals > 1
            no_intervals: it[0].num_intervals <= 1
        }.set{mpileup_normal}

    MPILEUP_TUMOR(cram_tumor, fasta)

    MPILEUP_TUMOR.out.mpileup.branch{
            intervals:    it[0].num_intervals > 1
            no_intervals: it[0].num_intervals <= 1
        }.set{mpileup_tumor}
        //TODO: Strelka anschauen? Single wird sowohl für tumor only und dingens verwendet
    //Merge mpileup only when intervals and natural order sort them
    CAT_MPILEUP_NORMAL(mpileup_normal.intervals
        .map{ meta, pileup ->
            new_meta = [patient:meta.patient, normal_id:meta.normal_id, tumor_id:meta.tumor_id, gender:meta.gender, id:meta.tumor_id + "_vs_" + meta.normal_id, num_intervals:meta.num_intervals]
            [groupKey(new_meta, meta.num_intervals), pileup]
            }
        .groupTuple(sort:true))

    CAT_MPILEUP_TUMOR(mpileup_tumor.intervals
        .map{ meta, pileup ->
            new_meta = [patient:meta.patient, normal_id:meta.normal_id, tumor_id:meta.tumor_id, gender:meta.gender, id:meta.tumor_id + "_vs_" + meta.normal_id, num_intervals:meta.num_intervals]
            [groupKey(new_meta, meta.num_intervals), pileup]
            }
        .groupTuple(sort:true))


    ch_versions = ch_versions.mix(MPILEUP_NORMAL.out.versions)
    ch_versions = ch_versions.mix(MPILEUP_TUMOR.out.versions)
    ch_versions = ch_versions.mix(CAT_MPILEUP_NORMAL.out.versions)
    ch_versions = ch_versions.mix(CAT_MPILEUP_TUMOR.out.versions)
    emit:
    versions = ch_versions
    cat_mpileup_normal = CAT_MPILEUP_NORMAL.out.file_out
    cat_mpileup_tumor = CAT_MPILEUP_TUMOR.out.file_out
    mpileup_normal_no_intervals = mpileup_normal.no_intervals
    mpileup_tumor_no_intervals = mpileup_tumor.no_intervals
}
