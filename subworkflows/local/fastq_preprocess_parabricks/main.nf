include { PARABRICKS_FQ2BAM               } from '../../../modules/nf-core/parabricks/fq2bam/main.nf'
include { CHANNEL_ALIGN_CREATE_CSV        } from '../../../subworkflows/local/channel_align_create_csv/main'
include { SAMTOOLS_CONVERT as CRAM_TO_BAM } from '../../../modules/nf-core/samtools/convert/main'
include { CRAM_MERGE_INDEX_SAMTOOLS       } from '../../../subworkflows/local/cram_merge_index_samtools/main'

workflow FASTQ_PREPROCESS_PARABRICKS {

    take:
    ch_reads                        // channel: [mandatory] meta, reads
    ch_fasta                        // channel: [mandatory] meta, fasta
    ch_fasta_fai                    // channel: [mandatory] meta, fasta_fai
    ch_index                        // channel: [mandatory] meta, index - bwa index
    ch_interval_file                // channel: [optional]  intervals_bed_combined
    ch_known_sites                  // channel: [optional]  known_sites_indels
    val_output_fmt                  // either bam or cram
    val_save_mapped                 // boolean
    val_save_output_as_bam          // boolean
    val_outdir                      // output directory for saving mapped files

    main:
    ch_versions = channel.empty()
    ch_reports  = channel.empty()

    ch_reads.map { meta, reads ->
            [ meta.subMap('patient', 'sample', 'sex', 'status'), reads ]
        }
        .groupTuple()
        .map { meta, reads ->
            meta + [ n_fastq: reads.size() ] // We can drop the FASTQ files now that we know how many there are
        }.set { reads_grouping_key }

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

    // Grouping the bams from the same samples not to stall the workflow
    // Use groupKey to make sure that the correct group can advance as soon as it is complete
    // and not stall the workflow until all reads from all channels are mapped
    cram_mapped = PARABRICKS_FQ2BAM.out.cram
        .combine(reads_grouping_key) // Creates a tuple of [ meta, bam, reads_grouping_key ]
        .filter { meta1, _cram, meta2 -> meta1.sample == meta2.sample }
        // Add n_fastq and other variables to meta
        .map { meta1, cram, meta2 ->
            [ meta1 + meta2, cram ]
        }
        // Manipulate meta map to remove old fields and add new ones
        .map { meta, cram ->
            [ meta - meta.subMap('id', 'read_group', 'data_type', 'num_lanes', 'size', 'sample_lane_id', 'lane') + [ data_type: 'cram', id: meta.sample ], cram ]
        }
        // Create groupKey from meta map
        .map { meta, cram ->
            [ groupKey( meta, meta.n_fastq), cram ]
        }
        // Group
        .groupTuple()

    // crams are merged (when multiple lanes from the same sample) and indexed
    CRAM_MERGE_INDEX_SAMTOOLS(cram_mapped, ch_fasta, ch_fasta_fai)

    ch_versions = ch_versions.mix(CRAM_MERGE_INDEX_SAMTOOLS.out.versions)

    cram_variant_calling = CRAM_MERGE_INDEX_SAMTOOLS.out.cram_crai
        .map { meta, cram, crai ->
                    [ meta - meta.subMap('id', 'read_group', 'data_type', 'num_lanes', 'size', 'sample_lane_id', 'lane') + [ data_type: 'cram', id: meta.sample ], cram, crai ]
                }

    if (val_save_output_as_bam) {
        // Convert CRAM files to BAM
        CRAM_TO_BAM(cram_variant_calling, ch_fasta, ch_fasta_fai)
        ch_versions = ch_versions.mix(CRAM_TO_BAM.out.versions)
        CHANNEL_ALIGN_CREATE_CSV(CRAM_TO_BAM.out.bam.join(CRAM_TO_BAM.out.bai, failOnDuplicate: true, failOnMismatch: true), val_outdir, val_save_output_as_bam)
    } else if (val_save_mapped) {
        CHANNEL_ALIGN_CREATE_CSV(cram_variant_calling, val_outdir, val_save_output_as_bam)
    }

    emit:
    cram      = cram_variant_calling     // channel: [ val(meta), cram, crai ]
    versions  = ch_versions              // channel: [ versions.yml ]
    reports   = ch_reports
}
