include { CAT_CAT as CAT              } from '../../../../modules/nf-core/modules/cat/cat/main.nf'
include { CONTROLFREEC                } from '../../../../modules/nf-core/modules/controlfreec/main'
include { SAMTOOLS_MPILEUP as MPILEUP } from '../../../../modules/nf-core/modules/samtools/mpileup/main'

workflow RUN_CONTROLFREEC {
    take:
    cram                     // channel: [mandatory] [meta, cram, crai, interval]
    fasta                    // channel: [mandatory]
    fasta_fai                // channel: [mandatory]
    dbsnp // channel: [mand]
    dbsnp_tbi
    chr_length
    mappability
    intervals_bed            // channel: [optional]  Contains a bed file of all intervals combined provided with the cram input(s). Should be empty for WGS
    num_intervals            //     val: [optional]  Number of used intervals, mandatory when intervals are provided.

    main:

    ch_versions = Channel.empty()

    MPILEUP(cram, fasta)

    MPILEUP.out.mpileup.branch{
            intervals:    num_intervals > 1
            no_intervals: num_intervals == 1
        }.set{mpileup}

    //Merge mpileup only when intervals
    CAT(mpileup.intervals.map{ meta, pileup ->
                new_meta = meta.clone()
                new_meta.id = new_meta.sample
                [new_meta, pileup]
            }.groupTuple(size: num_intervals))

    controlfeec_input = Channel.empty().mix(
        CAT.out.file_out,
        mpileup.no_intervals
    ).map{ meta, pileup ->
        [meta, pileup, [], [], [], [], []]
    }

    CONTROLFREEC(controlfeec_input,
                fasta,
                fasta_fai,
                [],
                dbsnp,
                dbsnp_tbi,
                chr_length,
                mappability,
                intervals_bed,
                [])
    //CONTROLFREECVis


    emit:
    versions = ch_versions
}
