//
// BAM TO CRAM and optionnal QC
//
// For all modules here:
// A when clause condition is defined in the conf/modules.config to determine if the module should be run

include { DEEPTOOLS_BAMCOVERAGE } from '../../modules/nf-core/modules/deeptools/bamcoverage/main'
include { QUALIMAP_BAMQC        } from '../../modules/nf-core/modules/qualimap/bamqc/main'
include { SAMTOOLS_BAMTOCRAM    } from '../../modules/nf-core/modules/samtools/bamtocram/main'

workflow BAM_TO_CRAM {
    take:
        bam_indexed                   // channel: [mandatory] meta, bam, bai
        fasta                         // channel: [mandatory] fasta
        fasta_fai                     // channel: [mandatory] fasta_fai
        intervals_combined_bed_gz_tbi // channel: [optional]  intervals_bed.gz, intervals_bed.gz.tbi

    main:
    ch_versions = Channel.empty()
    qc_reports  = Channel.empty()

    // remap to have channel without bam index
    bam_no_index = bam_indexed.map{ meta, bam, bai -> [meta, bam] }

    // Convert bam input to cram
    SAMTOOLS_BAMTOCRAM(bam_indexed, fasta)

    // Reports on bam input
    DEEPTOOLS_BAMCOVERAGE(bam_indexed)
    QUALIMAP_BAMQC(bam_no_index, intervals_combined_bed_gz_tbi)

    // Other reports run on cram

    // Gather all reports generated
    qc_reports = qc_reports.mix(DEEPTOOLS_BAMCOVERAGE.out.bigwig)
    qc_reports = qc_reports.mix(QUALIMAP_BAMQC.out.results)

    // Gather versions of all tools used
    ch_versions = ch_versions.mix(DEEPTOOLS_BAMCOVERAGE.out.versions.first())
    ch_versions = ch_versions.mix(QUALIMAP_BAMQC.out.versions.first())
    ch_versions = ch_versions.mix(SAMTOOLS_BAMTOCRAM.out.versions.first())

    emit:
        cram     = SAMTOOLS_BAMTOCRAM.out.cram_crai
        qc       = qc_reports

        versions = ch_versions // channel: [ versions.yml ]
}
