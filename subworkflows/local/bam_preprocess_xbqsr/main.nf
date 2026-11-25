include { SAMTOOLS_CONVERT as CRAM_TO_BAM } from '../../../modules/nf-core/samtools/convert/main'
include { SAMTOOLS_CONVERT as BAM_TO_CRAM } from '../../../modules/nf-core/samtools/convert/main'
include { SAMTOOLS_INDEX as BAM_INDEX     } from '../../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_INDEX as BAM_INDEX_RECAL } from '../../../modules/nf-core/samtools/index/main'
include { XBQSR_RECALIBRATE               } from '../../../modules/local/xbqsr/recalibrate/main'
include { XBQSR_APPLY                     } from '../../../modules/local/xbqsr/apply/main'

workflow BAM_PREPROCESS_XBQSR {
    take:
    input              // channel: [mandatory] [ meta, file, index ]
    dict               // channel: [mandatory] [ dict ]
    fasta              // channel: [mandatory] [ fasta ]
    fasta_fai          // channel: [mandatory] [ fasta_fai ]
    known_sites        // channel: [mandatory] [ known_sites ]
    known_sites_tbi    // channel: [optional]  [ known_sites_tbi ]
    save_output_as_bam // boolean: params.save_output_as_bam

    main:
    ch_versions = channel.empty()

    // Reference tuple for SAMTOOLS_CONVERT: [ meta, fasta, fai ]
    fasta_and_fai = fasta.combine(fasta_fai).map { meta_f, fasta_, _meta_fai, fai -> [ meta_f, fasta_, fai ] }.collect()

    // Branch inputs by file format (CRAM vs BAM)
    ch_input_branch = input.branch { meta, file, index ->
        cram: file.name.endsWith('.cram')
        bam:  file.name.endsWith('.bam')
    }

    // Convert CRAM to BAM if needed
    CRAM_TO_BAM(ch_input_branch.cram, fasta_and_fai)
    ch_versions = ch_versions.mix(CRAM_TO_BAM.out.versions_samtools)

    ch_converted_bam = CRAM_TO_BAM.out.bam.join(CRAM_TO_BAM.out.bai, failOnMismatch: true)

    // Ensure BAM has index
    ch_bam_to_index = ch_input_branch.bam.filter { meta, file, index -> !index }.map { meta, file, _idx -> [meta, file] }
    BAM_INDEX(ch_bam_to_index)
    ch_indexed_bam = ch_bam_to_index.join(BAM_INDEX.out.index, failOnMismatch: true)

    ch_bam_ready = ch_converted_bam
        .mix(ch_input_branch.bam.filter { meta, file, index -> index })
        .mix(ch_indexed_bam)

    // Execute X-BQSR Recalibrate
    XBQSR_RECALIBRATE(
        ch_bam_ready,
        fasta,
        known_sites.map { files -> [ [id: 'known_sites'], files ] }.collect()
    )
    ch_versions = ch_versions.mix(XBQSR_RECALIBRATE.out.versions_bqsr)
    ch_table = XBQSR_RECALIBRATE.out.table

    // Execute X-BQSR Apply
    ch_apply_input = ch_bam_ready.join(ch_table, failOnMismatch: true)
    XBQSR_APPLY(ch_apply_input)
    ch_versions = ch_versions.mix(XBQSR_APPLY.out.versions_bqsr)

    // Post-process according to save_output_as_bam
    if (save_output_as_bam) {
        BAM_INDEX_RECAL(XBQSR_APPLY.out.bam)
        ch_versions = ch_versions.mix(BAM_INDEX_RECAL.out.versions_samtools)
        recal_out = XBQSR_APPLY.out.bam.join(BAM_INDEX_RECAL.out.index, failOnMismatch: true)
    } else {
        BAM_TO_CRAM(XBQSR_APPLY.out.bam.map { meta, bam -> [meta, bam, []] }, fasta_and_fai)
        ch_versions = ch_versions.mix(BAM_TO_CRAM.out.versions_samtools)
        recal_out = BAM_TO_CRAM.out.cram.join(BAM_TO_CRAM.out.crai, failOnMismatch: true)
    }

    emit:
    recal_out = recal_out                                   // channel: [ meta, file, index ]
    table     = ch_table                                    // channel: [ meta, table ]
    reports   = ch_table.collect { _meta, table_f -> [ table_f ] }
    versions  = ch_versions
}
