//
// QC on CRAM
//
// For all modules here:
// A when clause condition is defined in the conf/modules.config to determine if the module should be run

include { SAMTOOLS_STATS } from '../../../modules/nf-core/samtools/stats/main'
include { MOSDEPTH       } from '../../../modules/nf-core/mosdepth/main'
include { IRONQC         } from '../../../modules/local/ironqc/main'

workflow CRAM_QC_MOSDEPTH_SAMTOOLS {
    take:
    cram      // channel: [mandatory] [ meta, cram, crai ]
    fasta     // channel: [mandatory] [ meta, fasta ]
    fasta_fai // channel: [mandatory] [ meta, fasta_fai ]
    intervals // channel: [optional]  [ meta, bed ]

    main:
    versions = Channel.empty()
    reports = Channel.empty()

    if (params.use_ironqc) {

        // IRONQC replaces both SAMTOOLS_STATS and MOSDEPTH in a single pass
        // Input: [meta, cram, crai, bed]
        cram_with_bed = cram.combine(intervals.map { meta, bed -> [bed ?: []] })

        IRONQC(
            cram_with_bed,
            fasta,
            fasta_fai
        )

        // Gather all reports generated (same outputs as original tools)
        reports = reports.mix(IRONQC.out.stats)
        reports = reports.mix(IRONQC.out.global_dist)
        reports = reports.mix(IRONQC.out.summary)

        // Gather versions
        versions = versions.mix(IRONQC.out.versions)

    } else {

        // Original tools: SAMTOOLS_STATS + MOSDEPTH
        SAMTOOLS_STATS(cram, fasta)

        MOSDEPTH(cram.combine(intervals.map { meta, bed -> [bed ?: []] }), fasta)

        // Gather all reports generated
        reports = reports.mix(SAMTOOLS_STATS.out.stats)
        reports = reports.mix(MOSDEPTH.out.global_txt)
        reports = reports.mix(MOSDEPTH.out.regions_txt)

        // Gather versions of all tools used
        versions = versions.mix(MOSDEPTH.out.versions)
        versions = versions.mix(SAMTOOLS_STATS.out.versions)
    }

    emit:
    reports
    versions // channel: [ versions.yml ]
}
