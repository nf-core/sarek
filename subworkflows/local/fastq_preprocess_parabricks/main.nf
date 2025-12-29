include { PARABRICKS_FQ2BAM               } from '../../../modules/nf-core/parabricks/fq2bam/main.nf'
include { CHANNEL_ALIGN_CREATE_CSV        } from '../../../subworkflows/local/channel_align_create_csv/main'
include { SAMTOOLS_CONVERT as CRAM_TO_BAM } from '../../../modules/nf-core/samtools/convert/main'

workflow FASTQ_PREPROCESS_PARABRICKS {

    take:
    ch_reads                        // channel: [mandatory] meta, reads
    ch_fasta                        // channel: [mandatory] meta, fasta
    ch_fasta_fai                    // channel: [mandatory] meta, fasta.fai
    ch_index                        // channel: [mandatory] meta, index - bwa index
    ch_interval_file                // channel: [optional]  intervals_bed_combined
    ch_known_sites                  // channel: [optional]  known_sites_indels
    val_output_fmt                  // either bam or cram

    main:
    ch_versions = channel.empty()
    ch_reports  = channel.empty()

    ch_reads.map { meta, reads ->
            [ meta.subMap('patient', 'sample', 'sex', 'status'), reads ]
        }
        .groupTuple()
        .map { meta, reads ->
            meta + [ n_fastq: reads.size() ] // We can drop the FASTQ files now that we know how many there are
        }

    ch_reads = ch_reads.map{ meta, reads ->
        // Update meta.id to meta.sample no multiple lanes or splitted fastqs
        if (meta.size * meta.num_lanes == 1) [ meta + [ id:meta.sample ], reads ]
        else [ meta, reads ]
    }

    // Adjust ch_interval_file
    ch_interval_file = ch_interval_file.collect().map { files ->
        [['id': 'intervals'], files]
    }

    // Adjust ch_known_sites
    ch_known_sites= ch_known_sites.collect().map { files ->
        [['id': 'known_sites'], files]
    }

    PARABRICKS_FQ2BAM(
        ch_reads,           // channel: [ val(meta), reads ]
        ch_fasta,           // channel: [ val(meta), fasta ]
        ch_index,           // channel: [ val(meta), index ]
        ch_interval_file,   // channel: [ val(meta), interval_file ]
        ch_known_sites,     // channel: [ val(meta), known_sites ]
        val_output_fmt      // either bam or cram
    )

    ch_versions = ch_versions.mix(PARABRICKS_FQ2BAM.out.versions)

    cram_variant_calling = PARABRICKS_FQ2BAM.out.cram
        .join(PARABRICKS_FQ2BAM.out.crai, failOnDuplicate: true, failOnMismatch: true)
        .map { meta, cram, crai ->
                    [ meta - meta.subMap('id', 'read_group', 'data_type', 'num_lanes', 'read_group', 'size', 'sample_lane_id') + [ data_type: 'cram', id: meta.sample ], cram, crai ]
                }

    // If params.save_output_as_bam, then convert CRAM files to BAM
    CRAM_TO_BAM(cram_variant_calling, ch_fasta, ch_fasta_fai)

    if (params.save_output_as_bam) CHANNEL_ALIGN_CREATE_CSV(CRAM_TO_BAM.out.bam.join(CRAM_TO_BAM.out.bai, failOnDuplicate: true, failOnMismatch: true), params.outdir, params.save_output_as_bam)
    else CHANNEL_ALIGN_CREATE_CSV(cram_variant_calling, params.outdir, params.save_output_as_bam)

    emit:
    cram      = cram_variant_calling     // channel: [ val(meta), cram, crai ]
    versions  = ch_versions              // channel: [ versions.yml ]
    reports   = ch_reports
}
