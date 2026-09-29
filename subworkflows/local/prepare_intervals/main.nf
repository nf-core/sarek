//
// PREPARE INTERVALS
//

// Initialize channels based on params or indices that were just built
// For all modules here:
// A when clause condition is defined in the conf/modules.config to determine if the module should be run

include { CREATE_INTERVALS_BED } from '../../../modules/local/create_intervals_bed'
include { GATK4_INTERVALLISTTOBED } from '../../../modules/nf-core/gatk4/intervallisttobed'
include { GAWK as BUILD_INTERVALS } from '../../../modules/nf-core/gawk'
include { HTSLIB_BGZIPTABIX as TABIX_BGZIPTABIX_INTERVAL_SPLIT } from '../../../modules/nf-core/htslib/bgziptabix'
include { HTSLIB_BGZIPTABIX as TABIX_BGZIPTABIX_INTERVAL_COMBINED } from '../../../modules/nf-core/htslib/bgziptabix'
include { taskOutputOrNull } from '../utils_nfcore_sarek_pipeline'
include { PreparedIntervals } from './types'

workflow PREPARE_INTERVALS {
    take:
    fasta_fai // mandatory [ fasta_fai ]
    intervals // [ params.intervals ]
    no_intervals // [ params.no_intervals ]
    nucleotides_per_second
    outdir
    step

    main:
    intervals_bed = channel.empty() // List of [ bed, num_intervals ], one for each region
    intervals_bed_gz_tbi = channel.empty() // List of [ bed.gz, bed,gz.tbi, num_intervals ], one for each region
    intervals_combined = channel.empty() // Single bed file containing all intervals

    if (no_intervals) {
        file("${outdir}/no_intervals.bed").text = ""
        file("${outdir}/no_intervals.bed.gz").text = "no_intervals\n"
        file("${outdir}/no_intervals.bed.gz.tbi").text = "no_intervals\n"

        intervals_bed = channel.fromPath(file("${outdir}/no_intervals.bed")).map { bed -> [bed, 0] }
        intervals_bed_gz_tbi = channel.fromPath(files("${outdir}/no_intervals.bed.{gz,gz.tbi}")).collect().map { files -> [files, 0] }
        intervals_combined = channel.fromPath(file("${outdir}/no_intervals.bed")).map { bed -> [[id: bed.simpleName], bed] }
    }
    else if (step != 'annotate' && step != 'controlfreec') {
        // If no interval/target file is provided, then generated intervals from FASTA file
        if (!intervals) {
            BUILD_INTERVALS(fasta_fai, [], [])

            intervals_combined = BUILD_INTERVALS.out.output

            CREATE_INTERVALS_BED(intervals_combined.map { _meta, path -> path }, nucleotides_per_second)

            intervals_bed = CREATE_INTERVALS_BED.out.bed
        }
        else {
            intervals_combined = channel.fromPath(file(intervals)).map { bed -> [[id: bed.baseName], bed] }
            CREATE_INTERVALS_BED(file(intervals), nucleotides_per_second)

            intervals_bed = CREATE_INTERVALS_BED.out.bed

            // If interval file is not provided as .bed, but e.g. as .interval_list then convert to BED format
            if (intervals.endsWith(".interval_list")) {
                GATK4_INTERVALLISTTOBED(intervals_combined)
                intervals_combined = GATK4_INTERVALLISTTOBED.out.bed
            }
        }

        // Now for the intervals.bed the following operations are done:
        // 1. Intervals file is split up into multiple bed files for scatter/gather
        // 2. Each bed file is indexed

        // 1. Intervals file is split up into multiple bed files for scatter/gather & grouping together small intervals
        intervals_bed = intervals_bed
            .flatten()
            .map { intervalFile ->
                def duration = 0.0
                intervalFile.eachLine { line ->
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
            // Adding number of intervals as elements
            .map { files -> [files, files.size()] }
            .transpose()

        // 2. Create bed.gz and bed.gz.tbi for each interval file. They are split by region (see above)
        TABIX_BGZIPTABIX_INTERVAL_SPLIT(intervals_bed.map { file, _num_intervals -> [[id: file.baseName], file, [], []] }, 'compress', true, 'bed')

        intervals_bed_gz_tbi = TABIX_BGZIPTABIX_INTERVAL_SPLIT.out.output
            .join(TABIX_BGZIPTABIX_INTERVAL_SPLIT.out.index)
            .map { _meta, bed, tbi -> [bed, tbi] }
            .toList()
            // Adding number of intervals as elements
            .map { files -> [files, files.size()] }
            .transpose()
    }

    TABIX_BGZIPTABIX_INTERVAL_COMBINED(intervals_combined.map { meta, bed -> [meta, bed, [], []] }, 'compress', true, 'bed')

    intervals_bed_combined = intervals_combined.map { _meta, bed -> bed }.collect()
    intervals_bed_gz_tbi_combined = TABIX_BGZIPTABIX_INTERVAL_COMBINED.out.output.join(TABIX_BGZIPTABIX_INTERVAL_COMBINED.out.index).map { _meta, gz, tbi -> [gz, tbi] }.collect()

    ch_results = intervals_bed
        .map { bed, _num_intervals -> bed }
        .toList()
        .map { items -> [items.collect { item -> taskOutputOrNull(item) }.findAll { item -> item != null } ?: null] }
        .combine(
            intervals_bed_gz_tbi.map { files, _num_intervals -> files[0] }.toList().map { items -> [items.collect { item -> taskOutputOrNull(item) }.findAll { item -> item != null } ?: null] }
        )
        .combine(intervals_bed_combined.map { items -> [taskOutputOrNull(items[0])] })
        .combine(intervals_bed_gz_tbi_combined.map { items -> [taskOutputOrNull(items[0])] })
        .filter { fields -> fields.any { field -> field != null } }
        .map { split_bed, split_bed_gz, combined_bed, combined_bed_gz ->
            record(
                split_bed: split_bed,
                split_bed_gz: split_bed_gz,
                combined_bed: combined_bed,
                combined_bed_gz: combined_bed_gz,
            )
        }

    emit:
    // Intervals split for parallel execution
    intervals_bed // [ intervals.bed, num_intervals ]
    intervals_bed_gz_tbi // [ intervals.bed.gz, intervals.bed.gz.tbi, num_intervals ]
    // All intervals in one file
    intervals_bed_combined // [ intervals.bed ]
    intervals_bed_gz_tbi_combined // [ intervals.bed.gz, intervals.bed.gz.tbi]
    results = ch_results // PreparedIntervals
}
