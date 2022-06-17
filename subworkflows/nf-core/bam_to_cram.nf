//
// BAM TO CRAM and optionnal QC
//
// For all modules here:
// A when clause condition is defined in the conf/modules.config to determine if the module should be run

include { DEEPTOOLS_BAMCOVERAGE                  } from '../../modules/nf-core/modules/deeptools/bamcoverage/main'
include { QUALIMAP_BAMQCCRAM                     } from '../../modules/nf-core/modules/qualimap/bamqccram/main'
include { SAMTOOLS_CONVERT as SAMTOOLS_BAMTOCRAM } from '../../modules/nf-core/modules/samtools/convert/main'
include { SAMTOOLS_STATS as SAMTOOLS_STATS_CRAM  } from '../../modules/nf-core/modules/samtools/stats/main'
include { MOSDEPTH } from '../../modules/nf-core/modules/mosdepth/main'

workflow BAM_TO_CRAM {
    take:
        bam_indexed                   // channel: [mandatory] meta, bam, bai
        cram_indexed
        fasta                         // channel: [mandatory] fasta
        fasta_fai                     // channel: [mandatory] fai
        intervals_combined_bed_gz_tbi // channel: [optional]  intervals_bed.gz, intervals_bed.gz.tbi
        intervals_bed_combined

    main:
    ch_versions = Channel.empty()
    qc_reports  = Channel.empty()

    // remap to have channel without bam index
    bam_no_index = bam_indexed.map{ meta, bam, bai -> [meta, bam] }

    // Convert bam input to cram
    SAMTOOLS_BAMTOCRAM(bam_indexed, fasta, fasta_fai)

    cram_indexed = Channel.empty().mix(cram_indexed,SAMTOOLS_BAMTOCRAM.out.alignment_index)

    // Reports on cram
    DEEPTOOLS_BAMCOVERAGE(cram_indexed, fasta, fasta_fai)
    QUALIMAP_BAMQCCRAM(cram_indexed, intervals_combined_bed_gz_tbi, fasta, fasta_fai)
    SAMTOOLS_STATS_CRAM(cram_indexed, fasta)
    MOSDEPTH(cram_indexed, intervals_bed_combined, fasta)

    // Gather all reports generated
    qc_reports = qc_reports.mix(DEEPTOOLS_BAMCOVERAGE.out.bigwig)
    qc_reports = qc_reports.mix(QUALIMAP_BAMQCCRAM.out.results)
    qc_reports = qc_reports.mix(SAMTOOLS_STATS_CRAM.out.stats)
    qc_reports = qc_reports.mix(MOSDEPTH.out.global_txt, MOSDEPTH.out.summary_txt, MOSDEPTH.out.per_base_d4)

    // Gather versions of all tools used
    ch_versions = ch_versions.mix(DEEPTOOLS_BAMCOVERAGE.out.versions.first())
    ch_versions = ch_versions.mix(QUALIMAP_BAMQCCRAM.out.versions.first())
    ch_versions = ch_versions.mix(SAMTOOLS_BAMTOCRAM.out.versions.first())
    ch_versions = ch_versions.mix(SAMTOOLS_STATS_CRAM.out.versions)

    emit:
        cram_converted  = SAMTOOLS_BAMTOCRAM.out.alignment_index
        qc              = qc_reports

        versions = ch_versions // channel: [ versions.yml ]
}
