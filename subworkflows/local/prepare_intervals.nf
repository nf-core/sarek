//
// PREPARE INTERVALS
//

// Initialize channels based on params or indices that were just built
// For all modules here:
// A when clause condition is defined in the conf/modules.config to determine if the module should be run

include { BUILD_INTERVALS                                     } from '../../modules/local/build_intervals/main'
include { CREATE_INTERVALS_BED                                } from '../../modules/local/create_intervals_bed/main'
include { GATK4_INTERVALLISTTOBED                             } from '../../modules/nf-core/modules/gatk4/intervallisttobed/main'
include { TABIX_BGZIPTABIX as TABIX_BGZIPTABIX_INTERVAL_SPLIT } from '../../modules/nf-core/modules/tabix/bgziptabix/main'
include { TABIX_BGZIPTABIX as TABIX_BGZIPTABIX_INTERVAL_ALL   } from '../../modules/nf-core/modules/tabix/bgziptabix/main'

workflow PREPARE_INTERVALS {
    take:
        fasta_fai // channel: [mandatory] fasta_fai

    main:

    ch_versions = Channel.empty()

    ch_intervals                     = Channel.empty()
    ch_intervals_bed_gz_tbi          = Channel.empty()
    ch_intervals_combined_bed_gz_tbi = Channel.empty() // Create bed.gz and bed.gz.tbi for input/or created interval file. Contains ALL regions.

    if (params.no_intervals) {
        file("${params.outdir}/no_intervals.bed").text = "no_intervals\n"
        file("${params.outdir}/no_intervals.bed.gz").text = "no_intervals\n"
        file("${params.outdir}/no_intervals.bed.gz.tbi").text = "no_intervals\n"

        ch_intervals = Channel.fromPath(file("${params.outdir}/no_intervals.bed"))
                                            .map{ it -> [it, 0]}

        ch_intervals_bed_gz_tbi = Channel.fromPath(file("${params.outdir}/no_intervals.bed.{gz,gz.tbi}"))
                                            .collect().map{ it -> [it, 0]}

        ch_intervals_combined_bed_gz_tbi = Channel.fromPath(file("${params.outdir}/no_intervals.bed.{gz,gz.tbi}"))
                                            .collect()

    } else if (params.step != 'annotate' && params.step != 'controlfreec') {

        tabix_in_combined                = Channel.empty()

        //If no interval/target file is provided, then intervals are generated from FASTA file
        if (!params.intervals) {

            BUILD_INTERVALS(fasta_fai)
            tabix_in_combined = BUILD_INTERVALS.out.bed.map{it -> [[id:it.simpleName], it] }

            ch_intervals = CREATE_INTERVALS_BED(BUILD_INTERVALS.out.bed)

        } else {

            tabix_in_combined = Channel.fromPath(file(params.intervals)).map{it -> [[id:it.baseName], it] }

            //If interval file is not provided as .bed, but e.g. as .interval_list then convert to BED format
            if(!params.intervals.endsWith(".bed")) {
                GATK4_INTERVALLISTTOBED(tabix_in_combined)
                tabix_in_combined = GATK4_INTERVALLISTTOBED.out.bed
                ch_versions = ch_versions.mix(GATK4_INTERVALLISTTOBED.out.versions)
            }

            ch_intervals = CREATE_INTERVALS_BED(file(params.intervals))
        }

        // Now for the interval.bed the following operations are done:
        // 1. Complete intervals file (with all intervals) is indexed
        // 2. Interval file is split up into multiple bed files for scatter/gather
        // 3. Each bed file from 2. is indexed

        // 1. Index complete interval file
        TABIX_BGZIPTABIX_INTERVAL_ALL(tabix_in_combined)
        ch_intervals_combined_bed_gz_tbi = TABIX_BGZIPTABIX_INTERVAL_ALL.out.gz_tbi.map{ meta, bed, tbi -> [bed, tbi] }
        ch_versions = ch_versions.mix(TABIX_BGZIPTABIX_INTERVAL_ALL.out.versions)

        // 2. Interval file is split up into multiple bed files for scatter/gather & grouping together small intervals
        ch_intervals = ch_intervals.flatten()
            .map{ intervalFile ->
                def duration = 0.0
                for (line in intervalFile.readLines()) {
                    final fields = line.split('\t')
                    if (fields.size() >= 5) duration += fields[4].toFloat()
                    else {
                        start = fields[1].toInteger()
                        end = fields[2].toInteger()
                        duration += (end - start) / params.nucleotides_per_second
                    }
                }
                [duration, intervalFile]
            }.toSortedList({ a, b -> b[0] <=> a[0] })
            .flatten().collate(2)
            .map{duration, intervalFile -> intervalFile}
            .collect().map{ it ->
                   [it, it.size() ] // Adding number of intervals as elements
                }.transpose()

        // 3. Create bed.gz and bed.gz.tbi for each interval file. They are split by region (see above)
        tabix_in = ch_intervals.map{ file, num_intervals -> [[id:file.baseName], file] }
        TABIX_BGZIPTABIX_INTERVAL_SPLIT(tabix_in)
        ch_intervals_bed_gz_tbi = TABIX_BGZIPTABIX_INTERVAL_SPLIT.out.gz_tbi.map{ meta, bed, tbi -> [bed, tbi ]}.toList().map{
                                        it ->
                                        [it, it.size()] // Adding number of intervals as elements
                                    }.transpose()
        ch_versions = ch_versions.mix(TABIX_BGZIPTABIX_INTERVAL_SPLIT.out.versions)

    }
    //ch_intervals.view()
    //ch_intervals_bed_gz_tbi.view()
    //ch_intervals_combined_bed_gz_tbi.view()

    emit:
        intervals_bed                    = ch_intervals                 // path: intervals.bed, num_intervals                        [intervals split for parallel execution]
        intervals_bed_gz_tbi             = ch_intervals_bed_gz_tbi          // path: target.bed.gz, target.bed.gz.tbi, num_intervals     [intervals split for parallel execution]
        intervals_combined_bed_gz_tbi    = ch_intervals_combined_bed_gz_tbi // path: interval.bed.gz, interval.bed.gz.tbi [all intervals in one file]
        versions                         = ch_versions                      // channel: [ versions.yml ]
}
