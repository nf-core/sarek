include { CAT_CAT as CAT_MPILEUP         } from '../../../../modules/nf-core/modules/cat/cat/main'
include { SAMTOOLS_MPILEUP               } from '../../../../modules/nf-core/modules/samtools/mpileup/main'

workflow RUN_MPILEUP {
    take:
        cram                    // channel: [mandatory] [meta, cram, interval]
        fasta                   // channel: [mandatory]

    main:

    ch_versions = Channel.empty()

    SAMTOOLS_MPILEUP(cram, fasta)
//    SAMTOOLS_MPILEUP.out.mpileup.view()
// [[patient:test1, sample:sample2, gender:XX, status:1, id:sample2_chr21_25689498-46709983, data_type:cram, num_intervals:2], /home/owacker/git/sarek/work/14/8c87bae2eb76aa70bce6637d850f22/sample2_chr21_25689498-46709983.mpileup]
    mpileup = SAMTOOLS_MPILEUP.out.mpileup.branch{
            intervals:    it[0].num_intervals > 1
            no_intervals: it[0].num_intervals <= 1
        }

//    mpileup.intervals.view()
//[[patient:test1, sample:sample2, gender:XX, status:1, id:sample2_chr21_25689498-46709983, data_type:cram, num_intervals:2], /home/owacker/git/sarek/work/c9/6804ec7cda0006fd651f401cedc983/sample2_chr21_25689498-46709983.mpileup]
    //Merge mpileup only when intervals and natural order sort them
    CAT_MPILEUP(mpileup.intervals
        .map{ meta, pileup ->
            new_meta = meta.tumor_id ? [patient:meta.patient, normal_id:meta.normal_id, tumor_id:meta.tumor_id, gender:meta.gender, id:meta.tumor_id + "_vs_" + meta.normal_id, num_intervals:meta.num_intervals] // not annotated, so no variantcaller necessary
                                        : [patient:meta.patient, sample:meta.sample, status:meta.status, gender:meta.gender, id:meta.id, num_intervals:meta.num_intervals] //Is id:sample or id:meta.id correct?
            [groupKey(new_meta, meta.num_intervals), pileup]
            }
                .groupTuple(sort:true))

    //TODO: THIS DOES NOT PRODUCE ANY OUTPUT!!!
    CAT_MPILEUP.out.file_out.view()

    Channel.empty().mix(CAT_MPILEUP.out.file_out, mpileup.no_intervals).view()

    ch_versions = ch_versions.mix(SAMTOOLS_MPILEUP.out.versions)
    ch_versions = ch_versions.mix(CAT_MPILEUP.out.versions)
    emit:
    versions = ch_versions
    mpileup = Channel.empty().mix(CAT_MPILEUP.out.file_out, mpileup.no_intervals)
}
