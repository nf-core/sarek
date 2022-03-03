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
        fasta_fai         // channel: [optional]  fasta_fai

    main:

    ch_versions = Channel.empty()

    ch_intervals                     = Channel.empty()
    ch_intervals_bed_gz_tbi          = Channel.empty()
    ch_intervals_combined_bed_gz_tbi = Channel.empty() // Create bed.gz and bed.gz.tbi for input/or created interval file. Contains ALL regions.

    tabix_in_combined = Channel.empty()

    if (params.no_intervals) {
        file("${params.outdir}/no_intervals.bed").text = "no_intervals\n"
        ch_intervals = Channel.fromPath(file("${params.outdir}/no_intervals.bed"))
        tabix_in_combined = ch_intervals.map{it -> [[id:it.getName()], it] }
    } else if (params.step != 'annotate' && params.step != 'controlfreec') {
        if (!params.intervals) {
            BUILD_INTERVALS(fasta_fai)
            tabix_in_combined = BUILD_INTERVALS.out.bed.map{it -> [[id:it.getName()], it] }
            ch_intervals = CREATE_INTERVALS_BED(BUILD_INTERVALS.out.bed)
        } else {
            tabix_in_combined = Channel.fromPath(file(params.intervals)).map{it -> [[id:it.baseName], it] }
            if(!params.intervals.endsWith(".bed")) {
                GATK4_INTERVALLISTTOBED(tabix_in_combined)
                tabix_in_combined = GATK4_INTERVALLISTTOBED.out.bed
                ch_versions = ch_versions.mix(GATK4_INTERVALLISTTOBED.out.versions)
            }
            ch_intervals = CREATE_INTERVALS_BED(file(params.intervals))
        }
    }
    if (params.step != 'annotate' && params.step != 'controlfreec') {

        TABIX_BGZIPTABIX_INTERVAL_ALL(tabix_in_combined)
        ch_intervals_combined_bed_gz_tbi = TABIX_BGZIPTABIX_INTERVAL_ALL.out.gz_tbi.map{ meta, bed, tbi -> [bed, tbi] }
        ch_versions = ch_versions.mix(TABIX_BGZIPTABIX_INTERVAL_ALL.out.versions)

        if (!params.no_intervals) {
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
        }

        // Create bed.gz and bed.gz.tbi for each interval file. They are split by region (see above)
        tabix_in = ch_intervals.map{it -> [[id:it.baseName], it] }
        TABIX_BGZIPTABIX_INTERVAL_SPLIT(tabix_in)
        ch_intervals_bed_gz_tbi = TABIX_BGZIPTABIX_INTERVAL_SPLIT.out.gz_tbi.map{ meta, bed, tbi -> [bed, tbi] }
        ch_versions = ch_versions.mix(TABIX_BGZIPTABIX_INTERVAL_SPLIT.out.versions)
    }

    emit:
        intervals_bed                    = ch_intervals                     // path: intervals.bed                        [intervals split for parallel execution]
        intervals_bed_gz_tbi             = ch_intervals_bed_gz_tbi          // path: target.bed.gz, target.bed.gz.tbi     [intervals split for parallel execution]
        intervals_combined_bed_gz_tbi    = ch_intervals_combined_bed_gz_tbi // path: interval.bed.gz, interval.bed.gz.tbi [all intervals in one file]

        versions                         = ch_versions                      // channel: [ versions.yml ]
}
