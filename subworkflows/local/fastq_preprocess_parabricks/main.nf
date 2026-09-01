include { PARABRICKS_APPLYBQSR      } from '../../../modules/nf-core/parabricks/applybqsr/main.nf'
include { PARABRICKS_FQ2BAM         } from '../../../modules/nf-core/parabricks/fq2bam/main.nf'
include { CHANNEL_ALIGN_CREATE_CSV  } from '../../../subworkflows/local/channel_align_create_csv/main'
include { BAM_MERGE_INDEX_SAMTOOLS  } from '../../../subworkflows/local/bam_merge_index_samtools/main'
include { CRAM_MERGE_INDEX_SAMTOOLS } from '../../../subworkflows/local/cram_merge_index_samtools/main'

workflow FASTQ_PREPROCESS_PARABRICKS {
    take:
    ch_reads // channel: [mandatory] meta, reads
    ch_fasta // channel: [mandatory] meta, fasta
    ch_fasta_fai // channel: [mandatory] meta, fasta_fai
    ch_index // channel: [mandatory] meta, index - bwa index
    ch_interval_file // channel: [optional]  intervals_bed_combined
    ch_known_sites // channel: [optional]  known_sites_indels
    val_skip_applybqsr // boolean
    val_save_mapped // boolean
    val_save_output_as_bam // boolean
    val_outdir // output directory for saving mapped files

    main:
    ch_reports = channel.empty()

    output_fmt = val_save_output_as_bam ? 'bam' : 'cram'

    reads_grouping_key = ch_reads
        .map { meta, reads ->
            [meta.subMap('patient', 'sample', 'sex', 'status'), reads]
        }
        .groupTuple()
        .map { meta, reads ->
            meta + [n_fastq: reads.size()]
        }

    ch_reads = ch_reads.map { meta, reads ->
        // Update meta.id to meta.sample no multiple lanes or splitted fastqs
        if (meta.size * meta.num_lanes == 1) {
            [meta + [id: meta.sample], reads]
        }
        else {
            [meta, reads]
        }
    }

    // Adjust ch_interval_file
    ch_interval_file = ch_interval_file
        .collect()
        .map { files ->
            [['id': 'intervals'], files]
        }

    // Adjust ch_known_sites
    ch_known_sites = ch_known_sites
        .collect()
        .map { files ->
            [['id': 'known_sites'], files]
        }

    PARABRICKS_FQ2BAM(
        ch_reads,
        ch_fasta,
        ch_index,
        ch_interval_file,
        ch_known_sites,
        output_fmt,
    )

    fq2bam_out_aln = val_save_output_as_bam ? PARABRICKS_FQ2BAM.out.bam : PARABRICKS_FQ2BAM.out.cram
    fq2bam_out_idx = val_save_output_as_bam ? PARABRICKS_FQ2BAM.out.bai : PARABRICKS_FQ2BAM.out.crai

    if (val_skip_applybqsr) {
        mapped_out_aln = fq2bam_out_aln
    }
    else {
        aln_idx = fq2bam_out_aln.join(fq2bam_out_idx)

        PARABRICKS_APPLYBQSR(
            aln_idx,
            PARABRICKS_FQ2BAM.out.bqsr_table,
            ch_interval_file,
            ch_fasta,
            output_fmt,
        )

        mapped_out_aln = val_save_output_as_bam ? PARABRICKS_APPLYBQSR.out.bam : PARABRICKS_APPLYBQSR.out.cram
    }

    if (val_save_output_as_bam) {
        // Grouping the CRAMs from the same samples not to stall the workflow
        // Use groupKey to make sure that the correct group can advance as soon as it is complete
        // and not stall the workflow until all reads from all channels are mapped
        bam_mapped = mapped_out_aln
            .combine(reads_grouping_key)
            .filter { meta1, _bam, meta2 -> meta1.sample == meta2.sample }
            .map { meta1, bam, meta2 ->
                [meta1 + meta2, bam]
            }
            .map { meta, bam ->
                [meta - meta.subMap('id', 'read_group', 'data_type', 'size', 'sample_lane_id', 'lane') + [data_type: 'bam', id: meta.sample], bam]
            }
            .map { meta, bam ->
                [groupKey(meta, meta.n_fastq), bam]
            }
            .groupTuple()

        // crams are merged (when multiple lanes from the same sample) and indexed
        BAM_MERGE_INDEX_SAMTOOLS(bam_mapped)

        bam_variant_calling = BAM_MERGE_INDEX_SAMTOOLS.out.bam_bai.map { meta, bam, bai ->
            [meta - meta.subMap('id', 'read_group', 'data_type', 'size', 'sample_lane_id', 'lane') + [data_type: 'bam', id: meta.sample], bam, bai]
        }

        CHANNEL_ALIGN_CREATE_CSV(bam_variant_calling, val_outdir, val_save_output_as_bam, true)

        // use created bam for downstream tools
        cram_variant_calling = bam_variant_calling
    }
    else {
        // Grouping the CRAMs from the same samples not to stall the workflow
        // Use groupKey to make sure that the correct group can advance as soon as it is complete
        // and not stall the workflow until all reads from all channels are mapped
        cram_mapped = mapped_out_aln
            .combine(reads_grouping_key)
            .filter { meta1, _cram, meta2 -> meta1.sample == meta2.sample }
            .map { meta1, cram, meta2 ->
                [meta1 + meta2, cram]
            }
            .map { meta, cram ->
                [meta - meta.subMap('id', 'read_group', 'data_type', 'size', 'sample_lane_id', 'lane') + [data_type: 'cram', id: meta.sample], cram]
            }
            .map { meta, cram ->
                [groupKey(meta, meta.n_fastq), cram]
            }
            .groupTuple()

        // crams are merged (when multiple lanes from the same sample) and indexed
        CRAM_MERGE_INDEX_SAMTOOLS(cram_mapped, ch_fasta, ch_fasta_fai)

        cram_variant_calling = CRAM_MERGE_INDEX_SAMTOOLS.out.cram_crai.map { meta, cram, crai ->
            [meta - meta.subMap('id', 'read_group', 'data_type', 'size', 'sample_lane_id', 'lane') + [data_type: 'cram', id: meta.sample], cram, crai]
        }

        if (val_save_mapped) {
            CHANNEL_ALIGN_CREATE_CSV(cram_variant_calling, val_outdir, val_save_mapped, true)
        }
    }

    emit:
    cram    = cram_variant_calling // channel: [ val(meta), cram, crai ]
    reports = ch_reports
}
