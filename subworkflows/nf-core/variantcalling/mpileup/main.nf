include { CAT_CAT as CAT_MPILEUP                    } from '../../../../modules/nf-core/modules/cat/cat/main'
include { SAMTOOLS_MPILEUP as MPILEUP               } from '../../../../modules/nf-core/modules/samtools/mpileup/main'

workflow RUN_MPILEUP {
    take:
        cram                    // channel: [mandatory] [meta, cram, crai, interval]
        fasta                   // channel: [mandatory]

    main:

    ch_versions = Channel.empty()

    MPILEUP(cram, fasta)

    MPILEUP.out.mpileup.branch{
            intervals:    it[0].num_intervals > 1
            no_intervals: it[0].num_intervals <= 1
        }
            .set{mpileup}

    //Merge mpileup only when intervals and natural order sort them
    CAT_MPILEUP(mpileup.intervals
        .map{ meta, pileup ->
            new_meta = meta.tumor_id ? [patient:meta.patient, normal_id:meta.normal_id, tumor_id:meta.tumor_id, gender:meta.gender, id:meta.tumor_id + "_vs_" + meta.normal_id, num_intervals:meta.num_intervals] // not annotated, so no variantcaller necessary
                                        : [patient:meta.patient, sample:meta.sample, status:meta.status, gender:meta.gender, id:new_id, num_intervals:meta.num_intervals]
            [groupKey(new_meta, meta.num_intervals), pileup]
            }
                .groupTuple(sort:true))

    ch_versions = ch_versions.mix(MPILEUP.out.versions)
    ch_versions = ch_versions.mix(CAT_MPILEUP.out.versions)
    emit:
    versions = ch_versions
    mpileup = Channel.empty().mix(CAT_MPILEUP.out.file_out, mpileup.no_intervals)
}
