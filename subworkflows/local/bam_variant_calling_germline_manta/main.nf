include { GATK4_MERGEVCFS as MERGE_MANTA_DIPLOID      } from '../../../modules/nf-core/gatk4/mergevcfs/main'
include { GATK4_MERGEVCFS as MERGE_MANTA_SMALL_INDELS } from '../../../modules/nf-core/gatk4/mergevcfs/main'
include { GATK4_MERGEVCFS as MERGE_MANTA_SV           } from '../../../modules/nf-core/gatk4/mergevcfs/main'
include { MANTA_GERMLINE                              } from '../../../modules/nf-core/manta/germline/main'

// Seems to be the consensus on upstream modules implementation too
workflow BAM_VARIANT_CALLING_GERMLINE_MANTA {
    take:
    cram                     // channel: [mandatory] [meta, cram, crai, interval.bed.gz, interval.bed.gz.tbi]
    dict                     // channel: [optional]
    fasta                    // channel: [mandatory]
    fasta_fai                // channel: [mandatory]

    main:
    versions = Channel.empty()

    MANTA_GERMLINE(cram, fasta, fasta_fai)

    // Figure out if using intervals or no_intervals
    small_indels_vcf = MANTA_GERMLINE.out.candidate_small_indels_vcf.branch{
        intervals:    it[0].num_intervals > 1
        no_intervals: it[0].num_intervals <= 1
    }

    sv_vcf = MANTA_GERMLINE.out.candidate_sv_vcf.branch{
        intervals:    it[0].num_intervals > 1
        no_intervals: it[0].num_intervals <= 1
    }

    diploid_sv_vcf = MANTA_GERMLINE.out.diploid_sv_vcf.branch{
        intervals:    it[0].num_intervals > 1
        no_intervals: it[0].num_intervals <= 1
    }

    // Only when using intervals
    MERGE_MANTA_SMALL_INDELS(
        small_indels_vcf.intervals.map{ meta, vcf ->
            [ groupKey(
                meta.subMap('num_intervals', 'patient', 'sample', 'sex', 'status')
                    + [ id:meta.sample ],
                        meta.num_intervals),
            vcf ]
        }.groupTuple(),
        dict.map{ it -> [ [ id:it[0].baseName ], it ] })

    MERGE_MANTA_SV(
        sv_vcf.intervals.map{ meta, vcf ->
            [ groupKey(
                meta.subMap('num_intervals', 'patient', 'sample', 'sex', 'status')
                    + [ id:meta.sample ],
                        meta.num_intervals),
            vcf ]
        }.groupTuple(),
        dict.map{ it -> [ [ id:it[0].baseName ], it ] })

    MERGE_MANTA_DIPLOID(
        diploid_sv_vcf.intervals.map{ meta, vcf ->
            [ groupKey(
                meta.subMap('num_intervals', 'patient', 'sample', 'sex', 'status')
                    + [ id:meta.sample ],
                        meta.num_intervals),
            vcf ]
        }.groupTuple(),
        dict.map{ it -> [ [ id:it[0].baseName ], it ] })

    // Mix output channels for "no intervals" and "with intervals" results
    // Only diploid SV should get annotated
    vcf = Channel.empty().mix(MERGE_MANTA_DIPLOID.out.vcf, diploid_sv_vcf.no_intervals).map{ meta, vcf ->
        [ meta.subMap('num_intervals', 'patient', 'sample', 'sex', 'status')
            + [ id:meta.sample, variantcaller: 'manta' ],
            vcf ]
    }

    versions = versions.mix(MERGE_MANTA_DIPLOID.out.versions)
    versions = versions.mix(MERGE_MANTA_SMALL_INDELS.out.versions)
    versions = versions.mix(MERGE_MANTA_SV.out.versions)
    versions = versions.mix(MANTA_GERMLINE.out.versions)

    emit:
    vcf

    versions
}
