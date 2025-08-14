//
// PREPARE RECALIBRATION SPARK
//
// For all modules here:
// A when clause condition is defined in the conf/modules.config to determine if the module should be run

include { GATK4SPARK_BASERECALIBRATOR } from '../../../modules/nf-core/gatk4spark/baserecalibrator'
include { GATK4_GATHERBQSRREPORTS     } from '../../../modules/nf-core/gatk4/gatherbqsrreports'

workflow BAM_BASERECALIBRATOR_SPARK {
    take:
    cram            // channel: [mandatory] [ meta, cram_markduplicates, crai ]
    dict            // channel: [mandatory] [ dict ]
    fasta           // channel: [mandatory] [ fasta ]
    fasta_fai       // channel: [mandatory] [ fasta_fai ]
    intervals       // channel: [mandatory] [ intervals, num_intervals ] (or [ [], 0 ] if no intervals)
    known_sites     // channel: [optional]  [ known_sites ]
    known_sites_tbi // channel: [optional]  [ known_sites_tbi ]

    main:
    versions = Channel.empty()

    // Combine cram and intervals for spread and gather strategy
    // Move num_intervals to meta map
    cram_intervals = cram
        .combine(intervals)
        .map { meta, cram_, crai, intervals_, num_intervals -> [meta + [num_intervals: num_intervals], cram_, crai, intervals_] }

    // RUN BASERECALIBRATOR SPARK
    GATK4SPARK_BASERECALIBRATOR(
        cram_intervals,
        fasta.map { _meta, fasta_ -> [fasta_] },
        fasta_fai.map { _meta, fasta_fai_ -> [fasta_fai_] },
        dict.map { _meta, dict_ -> [dict_] },
        known_sites,
        known_sites_tbi,
    )

    // Figuring out if there is one or more table(s) from the same sample
    table_to_merge = GATK4SPARK_BASERECALIBRATOR.out.table
        .map { meta, table -> [groupKey(meta, meta.num_intervals), table] }
        .groupTuple()
        .branch {
            single: it[0].num_intervals <= 1
            multiple: it[0].num_intervals > 1
        }

    // Only when using intervals
    GATK4_GATHERBQSRREPORTS(table_to_merge.multiple)

    // Mix intervals and no_intervals channels together
    // Remove no longer necessary field: num_intervals
    table_bqsr = GATK4_GATHERBQSRREPORTS.out.table
        .mix(table_to_merge.single.map { meta, table -> [meta, table[0]] })
        .map { meta, table -> [meta - meta.subMap('num_intervals'), table] }

    // Gather versions of all tools used
    versions = versions.mix(GATK4SPARK_BASERECALIBRATOR.out.versions)
    versions = versions.mix(GATK4_GATHERBQSRREPORTS.out.versions)

    emit:
    table_bqsr // channel: [ meta, table ]
    versions   // channel: [ versions.yml ]
}
