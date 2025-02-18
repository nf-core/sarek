//
// PREPARE INTERVALS
//

// Initialize channels based on params or indices that were just built
// For all modules here:
// A when clause condition is defined in the conf/modules.config to determine if the module should be run

include { CREATE_INTERVALS_BED                                   } from '../../../modules/local/create_intervals_bed'
include { GATK4_INTERVALLISTTOBED                                } from '../../../modules/nf-core/gatk4/intervallisttobed'
include { TABIX_BGZIPTABIX as TABIX_BGZIPTABIX_INTERVAL_SPLIT    } from '../../../modules/nf-core/tabix/bgziptabix'
include { TABIX_BGZIPTABIX as TABIX_BGZIPTABIX_INTERVAL_COMBINED } from '../../../modules/nf-core/tabix/bgziptabix'

workflow PREPARE_INTERVALS {
    take:
    intervals              // mandatory [ intervals_bed ]
    no_intervals           // boolean   [ params.no_intervals ]
    nucleotides_per_second // mandatory [ params.nucleotides_per_second ]
    outdir                 // mandatory [ params.outdir ]
    step                   // mandatory [ params.step ]

    main:
    versions = Channel.empty()

    // intervals_bed          - List of [ bed, num_intervals ], one per region
    // intervals_bed_gz_tbi   - List of [ bed.gz, bed,gz.tbi, num_intervals ], one per region
    // intervals_bed_combined - Single bed file containing all intervals
    intervals_bed = Channel.empty()
    intervals_bed_gz_tbi = Channel.empty()
    intervals_bed_combined = Channel.empty()

    if (no_intervals) {
        file("${outdir}/no_intervals.bed").text = "no_intervals\n"
        file("${outdir}/no_intervals.bed.gz").text = "no_intervals\n"
        file("${outdir}/no_intervals.bed.gz.tbi").text = "no_intervals\n"

        intervals_bed = Channel.fromPath(file("${outdir}/no_intervals.bed")).map { it -> [it, 0] }
        intervals_bed_gz_tbi = Channel.fromPath(file("${outdir}/no_intervals.bed.{gz,gz.tbi}")).collect().map { it -> [it, 0] }
        intervals_bed_combined = Channel.fromPath(file("${outdir}/no_intervals.bed")).map { it -> [[id: it.simpleName], it] }
    }
    else if (step != 'annotate' && step != 'controlfreec') {
        CREATE_INTERVALS_BED(intervals, nucleotides_per_second)

        intervals_bed = CREATE_INTERVALS_BED.out.bed

        versions = versions.mix(CREATE_INTERVALS_BED.out.versions)

        intervals_branch = intervals.branch { _meta, intervals_ ->
            interval_list: intervals_.endsWith(".interval_list")
            bed: true
        }

        GATK4_INTERVALLISTTOBED(intervals_branch.interval_list)
        // TODO: test this with an interval_list
        intervals_bed_combined = intervals.mix(GATK4_INTERVALLISTTOBED.out.bed).last()
        versions = versions.mix(GATK4_INTERVALLISTTOBED.out.versions)

        // Now for the intervals.bed the following operations are done:
        // 1/ Split up intervals bed file into multiple bed files for scatter/gather
        // 2/ Tabix index each bed file

        // 1/ Split up intervals bed file into multiple bed files for scatter/gather
        //      Also group together small intervals
        //      And add the number of intervals to the channel
        intervals_bed = intervals_bed
            .flatten()
            .map { intervals_ ->
                def duration = 0.0
                intervals_
                    .readLines()
                    .each { line ->
                        def fields = line.split('\t')
                        if (fields.size() >= 5) {
                            duration += fields[4].toFloat()
                        }
                        else {
                            def start = fields[1].toInteger()
                            def end = fields[2].toInteger()
                            duration += (end - start) / nucleotides_per_second
                        }
                    }
                [duration, intervals_]
            }
            .toSortedList { a, b -> b[0] <=> a[0] }
            .flatten()
            .collate(2)
            .map { _duration, intervals_ -> intervals_ }
            .collect()
            .map { intervals_ -> [intervals_, intervals_.size()] }
            .transpose()

        // 2/ Tabix index each bed file
        TABIX_BGZIPTABIX_INTERVAL_SPLIT(intervals_bed.map { intervals_, _num_intervals -> [[id: intervals_.baseName], intervals_] })

        intervals_bed_gz_tbi = TABIX_BGZIPTABIX_INTERVAL_SPLIT.out.gz_tbi
            .map { _meta, intervals_gz_, intervals_gz_tbi_ -> [intervals_gz_, intervals_gz_tbi_] }
            .toList()
            .map { it -> [it, it.size()] }
            .transpose()

        versions = versions.mix(TABIX_BGZIPTABIX_INTERVAL_SPLIT.out.versions)
    }

    TABIX_BGZIPTABIX_INTERVAL_COMBINED(intervals_bed_combined)
    versions = versions.mix(TABIX_BGZIPTABIX_INTERVAL_COMBINED.out.versions)

    // intervals_bed and intervals_bed_gz_tbi are the intervals split for parallel execution, and contain the number of intervals
    // intervals_bed_combined and intervals_bed_gz_tbi_combined are all intervals collected in one file

    intervals_bed_combined = intervals_bed_combined.map { _meta, intervals_ -> intervals_ }.collect()
    intervals_bed_gz_tbi_combined = TABIX_BGZIPTABIX_INTERVAL_COMBINED.out.gz_tbi.map { _meta, intervals_gz, intervals_gz_tbi -> [intervals_gz, intervals_gz_tbi] }.collect()

    emit:
    intervals_bed                 // [ intervals.bed, num_intervals ]
    intervals_bed_gz_tbi          // [ intervals.bed.gz, intervals.bed.gz.tbi, num_intervals ]
    intervals_bed_combined        // [ intervals.bed ]
    intervals_bed_gz_tbi_combined // [ intervals.bed.gz, intervals.bed.gz.tbi]
    versions                      // [ versions.yml ]
}
