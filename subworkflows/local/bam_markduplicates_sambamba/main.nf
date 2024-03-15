//
// MARKDUPLICATES AND QC after mapping
//
// For all modules here:
// A when clause condition is defined in the conf/modules.config to determine if the module should be run

include { CRAM_QC_MOSDEPTH_SAMTOOLS } from '../cram_qc_mosdepth_samtools/main'
include { SAMBAMBA_MARKDUP      } from '../../../modules/nf-core/sambamba/markdup/main'

workflow BAM_MARKDUPLICATES_SAMBAMBA {
    take:
    bam                    // channel: [mandatory] [ meta, bam ]
    fasta                  // channel: [mandatory] [ fasta ]
    fasta_fai              // channel: [mandatory] [ fasta_fai ]
    intervals_bed_combined // channel: [optional]  [ intervals_bed ]

    main:
    versions = Channel.empty()
    reports  = Channel.empty()

    // RUN MARKUPDUPLICATES
    bam.view()
    SAMBAMBA_MARKDUP(bam)

    // bam = SAMBAMBA_MARKDUP.out.bam
    // bam.view()

    cram = bam
    // Join with the crai file
    // bam = SAMBAMBA_MARKDUP.out.bam.join(SAMBAMBA_MARKDUPLICATES.out.bai, failOnDuplicate: true, failOnMismatch: true)

    // QC on CRAM
    // CRAM_QC_MOSDEPTH_SAMTOOLS(bam, fasta, intervals_bed_combined)

    // Gather all reports generated
    // reports = reports.mix(SAMBAMBA_MARKDUP.out.metrics)
    // reports = reports.mix(CRAM_QC_MOSDEPTH_SAMTOOLS.out.reports)

    // Gather versions of all tools used
    versions = versions.mix(SAMBAMBA_MARKDUP.out.versions)
    // versions = versions.mix(CRAM_QC_MOSDEPTH_SAMTOOLS.out.versions)

    emit:
    cram
    reports

    versions    // channel: [ versions.yml ]
}
