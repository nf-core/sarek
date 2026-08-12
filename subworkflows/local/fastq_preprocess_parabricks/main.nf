include { PARABRICKS_APPLYBQSR                    } from '../../../modules/nf-core/parabricks/applybqsr/main.nf'
include { PARABRICKS_FQ2BAM                       } from '../../../modules/nf-core/parabricks/fq2bam/main.nf'
include { CHANNEL_ALIGN_CREATE_CSV                } from '../../../subworkflows/local/channel_align_create_csv/main'
include { SAMTOOLS_CONVERT as CRAM_TO_BAM         } from '../../../modules/nf-core/samtools/convert/main'
include { SAMTOOLS_CONVERT as CRAM_TO_BAM_TO_SAVE } from '../../../modules/nf-core/samtools/convert/main'
include { CRAM_MERGE_INDEX_SAMTOOLS               } from '../../../subworkflows/local/cram_merge_index_samtools/main'

workflow FASTQ_PREPROCESS_PARABRICKS {

    take:
    ch_reads                        // channel: [mandatory] meta, reads
    ch_fasta                        // channel: [mandatory] meta, fasta
    ch_fasta_fai                    // channel: [mandatory] meta, fasta_fai
    ch_index                        // channel: [mandatory] meta, index - bwa index
    ch_interval_file                // channel: [optional]  intervals_bed_combined
    ch_known_sites                  // channel: [optional]  known_sites_indels
    val_output_fmt                  // either bam or cram
    val_applybqsr                   // boolean
    val_save_mapped                 // boolean
    val_save_output_as_bam          // boolean
    val_outdir                      // output directory for saving mapped files

    main:
    ch_reports  = channel.empty()

    reads_grouping_key = ch_reads.map { meta, reads ->
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

    // Single reference genome used throughout the pipeline
    fasta_fai = ch_fasta.combine(ch_fasta_fai).map { meta, fasta_, _meta_fai, fai -> [ meta, fasta_, fai ] }.collect()

    // BQSR needs fq2bam to emit BAM (applybqsr only ingests BAM), otherwise we keep the requested format
    fq2bam_output_fmt = val_applybqsr? 'bam' : val_output_fmt

    PARABRICKS_FQ2BAM(
        ch_reads,           // channel: [ val(meta), reads ]
        ch_fasta,           // channel: [ val(meta), fasta ]
        ch_index,           // channel: [ val(meta), index ]
        ch_interval_file,   // channel: [ val(meta), interval_file ]
        ch_known_sites,     // channel: [ val(meta), known_sites ]
        fq2bam_output_fmt   // either bam or cram
    )

    if (val_applybqsr) {
        // Combine each fq2bam BAM with its index, recalibration table (same run -> same
        // meta), the optional intervals and the reference into a single tuple channel.
        // The applybqsr module takes everything in ONE tuple so this works for any number
        // of samples (a plain multi-tuple process input would not broadcast the singletons).
        bqsr_input = PARABRICKS_FQ2BAM.out.bam
            .join(PARABRICKS_FQ2BAM.out.bai, failOnDuplicate: true, failOnMismatch: true)
            .join(PARABRICKS_FQ2BAM.out.bqsr_table, failOnDuplicate: true, failOnMismatch: true)
            .combine(ch_interval_file)
            .combine(ch_fasta)
            .map { meta, bam_, bai, table, _imeta, intervals, _fmeta, fasta -> [ meta, bam_, bai, table, intervals, fasta ] }

        PARABRICKS_APPLYBQSR(bqsr_input)

        // Native BAM (per lane) is the applybqsr output to keep when saving as BAM
        native_bam = PARABRICKS_APPLYBQSR.out.bam.join(PARABRICKS_APPLYBQSR.out.bai, failOnDuplicate: true, failOnMismatch: true)

        // Inverted vs the non-BQSR path (where fq2bam already yields CRAM): the recalibrated
        // BAM is converted back to CRAM so the rest of the pipeline still handles CRAMs.
        CRAM_TO_BAM(native_bam, fasta_fai)

        mapped_cram = CRAM_TO_BAM.out.cram
    }
    else {
        mapped_cram = PARABRICKS_FQ2BAM.out.cram
        native_bam  = channel.empty()
    }

    // Grouping the bams from the same samples not to stall the workflow
    // Use groupKey to make sure that the correct group can advance as soon as it is complete
    // and not stall the workflow until all reads from all channels are mapped
    cram_mapped = mapped_cram
        .combine(reads_grouping_key) // Creates a tuple of [ meta, bam, reads_grouping_key ]
        .filter { meta1, _cram, meta2 -> meta1.sample == meta2.sample }
        // Add n_fastq and other variables to meta
        .map { meta1, cram, meta2 ->
            [ meta1 + meta2, cram ]
        }
        // Manipulate meta map to remove old fields and add new ones
        .map { meta, cram ->
                [ meta - meta.subMap('id', 'read_group', 'data_type', 'size', 'sample_lane_id', 'lane') + [ data_type: 'cram', id: meta.sample ], cram ]
        }
        // Create groupKey from meta map
        .map { meta, cram ->
            [ groupKey( meta, meta.n_fastq), cram ]
        }
        // Group
        .groupTuple()

    // crams are merged (when multiple lanes from the same sample) and indexed
    CRAM_MERGE_INDEX_SAMTOOLS(cram_mapped, ch_fasta, ch_fasta_fai)

    cram_variant_calling = CRAM_MERGE_INDEX_SAMTOOLS.out.cram_crai
        .map { meta, cram, crai ->
                [ meta - meta.subMap('id', 'read_group', 'data_type', 'size', 'sample_lane_id', 'lane') + [ data_type: 'cram', id: meta.sample ], cram, crai ]
            }

    if (val_save_output_as_bam) {
        if (val_applybqsr) {
            // BQSR already produced BAMs, so saving as BAM needs no conversion (inverted logic)
            CHANNEL_ALIGN_CREATE_CSV(native_bam, val_outdir, val_save_output_as_bam)
        } else {
            // Convert CRAM files to BAM
            CRAM_TO_BAM_TO_SAVE(cram_variant_calling, fasta_fai)
            CHANNEL_ALIGN_CREATE_CSV(CRAM_TO_BAM_TO_SAVE.out.bam.join(CRAM_TO_BAM_TO_SAVE.out.bai, failOnDuplicate: true, failOnMismatch: true), val_outdir, val_save_output_as_bam)
        }
    } else if (val_save_mapped) {
        CHANNEL_ALIGN_CREATE_CSV(cram_variant_calling, val_outdir, val_save_output_as_bam)
    }

    emit:
    cram      = cram_variant_calling     // channel: [ val(meta), cram, crai ]
    reports   = ch_reports
}
