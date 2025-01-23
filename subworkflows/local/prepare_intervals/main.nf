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
    intervals_bed          // [ meta, intervals_bed ]
    no_intervals           // [ params.no_intervals ]
    nucleotides_per_second // [ params.nucleotides_per_second ]
    outdir                 // [ params.outdir ]
    step                   // [ params.step ]

    main:
    versions = Channel.empty()

    // intervals_bed        - List of [ bed, num_intervals ], one for each region
    // intervals_bed_gz_tbi - List of [ bed.gz, bed,gz.tbi, num_intervals ], one for each region
    // intervals_combined   - Single bed file containing all intervals
    intervals_bed = Channel.empty()
    intervals_bed_gz_tbi = Channel.empty()
    intervals_combined = Channel.empty()

    if (no_intervals) {
        file("${outdir}/no_intervals.bed").text = "no_intervals\n"
        file("${outdir}/no_intervals.bed.gz").text = "no_intervals\n"
        file("${outdir}/no_intervals.bed.gz.tbi").text = "no_intervals\n"

        intervals_bed = Channel.fromPath(file("${outdir}/no_intervals.bed")).map { it -> [it, 0] }
        intervals_bed_gz_tbi = Channel.fromPath(file("${outdir}/no_intervals.bed.{gz,gz.tbi}")).collect().map { it -> [it, 0] }
        intervals_combined = Channel.fromPath(file("${outdir}/no_intervals.bed")).map { it -> [[id: it.simpleName], it] }
    }
    else if (step != 'annotate' && step != 'controlfreec') {
        intervals_combined = intervals_bed
        CREATE_INTERVALS_BED(intervals_bed, nucleotides_per_second)

        intervals_bed = CREATE_INTERVALS_BED.out.bed

        versions = versions.mix(CREATE_INTERVALS_BED.out.versions)

        // // If interval file is not provided as .bed, but e.g. as .interval_list then convert to BED format
        // if (intervals.endsWith(".interval_list")) {
        //     GATK4_INTERVALLISTTOBED(intervals_combined)
        //     intervals_combined = GATK4_INTERVALLISTTOBED.out.bed
        //     versions = versions.mix(GATK4_INTERVALLISTTOBED.out.versions)
        // }

        // Now for the intervals.bed the following operations are done:
        // 1. Intervals file is split up into multiple bed files for scatter/gather
        // 2. Each bed file is indexed

        // 1. Intervals file is split up into multiple bed files for scatter/gather & grouping together small intervals
        intervals_bed = intervals_bed
            .flatten()
            .map { intervalFile ->
                def duration = 0.0
                intervalFile
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
                [duration, intervalFile]
            }
            .toSortedList { a, b -> b[0] <=> a[0] }
            .flatten()
            .collate(2)
            .map { _duration, intervalFile -> intervalFile }
            .collect()
            .map { it -> [it, it.size()] }
            .transpose()

        // 2. Create bed.gz and bed.gz.tbi for each interval file. They are split by region (see above)
        TABIX_BGZIPTABIX_INTERVAL_SPLIT(intervals_bed.map { file, _num_intervals -> [[id: file.baseName], file] })

        intervals_bed_gz_tbi = TABIX_BGZIPTABIX_INTERVAL_SPLIT.out.gz_tbi
            .map { _meta, bed, tbi -> [bed, tbi] }
            .toList()
            .map { it -> [it, it.size()] }
            .transpose()

        versions = versions.mix(TABIX_BGZIPTABIX_INTERVAL_SPLIT.out.versions)
    }

    TABIX_BGZIPTABIX_INTERVAL_COMBINED(intervals_combined)
    versions = versions.mix(TABIX_BGZIPTABIX_INTERVAL_COMBINED.out.versions)

    // intervals_bed and intervals_bed_gz_tbi are the intervals split for parallel execution, and contain the number of intervals
    // intervals_bed_combined and intervals_bed_gz_tbi_combined are all intervals collected in one file

    intervals_bed_combined = intervals_combined.map { _meta, bed -> bed }.collect()
    intervals_bed_gz_tbi_combined = TABIX_BGZIPTABIX_INTERVAL_COMBINED.out.gz_tbi.map { _meta, gz, tbi -> [gz, tbi] }.collect()

    emit:
    intervals_bed                 // [ intervals.bed, num_intervals ]
    intervals_bed_gz_tbi          // [ intervals.bed.gz, intervals.bed.gz.tbi, num_intervals ]
    intervals_bed_combined        // [ intervals.bed ]
    intervals_bed_gz_tbi_combined // [ intervals.bed.gz, intervals.bed.gz.tbi]
    versions                      // [ versions.yml ]
}
