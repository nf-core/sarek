// Create samplesheets to restart from different steps
include { CHANNEL_ALIGN_CREATE_CSV                          } from '../../../subworkflows/local/channel_align_create_csv/main'
include { CHANNEL_MARKDUPLICATES_CREATE_CSV                 } from '../../../subworkflows/local/channel_markduplicates_create_csv/main'
include { CHANNEL_BASERECALIBRATOR_CREATE_CSV               } from '../../../subworkflows/local/channel_baserecalibrator_create_csv/main'
include { CHANNEL_APPLYBQSR_CREATE_CSV                      } from '../../../subworkflows/local/channel_applybqsr_create_csv/main'

// Convert BAM files to FASTQ files
include { BAM_CONVERT_SAMTOOLS as CONVERT_FASTQ_UMI         } from '../../../subworkflows/local/bam_convert_samtools/main'

// TRIM/SPLIT FASTQ Files
include { FASTP                                             } from '../../../modules/nf-core/fastp/main'

// Create umi consensus bams from fastq
include { FASTQ_CREATE_UMI_CONSENSUS_FGBIO                  } from '../../../subworkflows/local/fastq_create_umi_consensus_fgbio/main'

// Map input reads to reference genome
include { FASTQ_ALIGN                                       } from '../../../subworkflows/local/fastq_align/main'

// Merge and index BAM files (optional)
include { BAM_MERGE_INDEX_SAMTOOLS                          } from '../../../subworkflows/local/bam_merge_index_samtools/main'

// Convert BAM files
include { SAMTOOLS_CONVERT as BAM_TO_CRAM_MAPPING           } from '../../../modules/nf-core/samtools/convert/main'

// Convert CRAM files (optional)
include { SAMTOOLS_CONVERT as CRAM_TO_BAM                   } from '../../../modules/nf-core/samtools/convert/main'
include { SAMTOOLS_CONVERT as CRAM_TO_BAM_RECAL             } from '../../../modules/nf-core/samtools/convert/main'

// Mark Duplicates (+QC)
include { BAM_MARKDUPLICATES                                } from '../../../subworkflows/local/bam_markduplicates/main'
include { BAM_MARKDUPLICATES_SPARK                          } from '../../../subworkflows/local/bam_markduplicates_spark/main'
include { BAM_SENTIEON_DEDUP                                } from '../../../subworkflows/local/bam_sentieon_dedup/main'

// QC on CRAM
include { CRAM_QC_MOSDEPTH_SAMTOOLS as CRAM_QC_NO_MD        } from '../../../subworkflows/local/cram_qc_mosdepth_samtools/main'

// Create recalibration tables
include { BAM_BASERECALIBRATOR                              } from '../../../subworkflows/local/bam_baserecalibrator/main'
include { BAM_BASERECALIBRATOR_SPARK                        } from '../../../subworkflows/local/bam_baserecalibrator_spark/main'

// Create recalibrated cram files to use for variant calling (+QC)
include { BAM_APPLYBQSR                                     } from '../../../subworkflows/local/bam_applybqsr/main'
include { BAM_APPLYBQSR_SPARK                               } from '../../../subworkflows/local/bam_applybqsr_spark/main'

workflow FASTQ_PREPROCESS_GATK {
    take:
        input_fastq
        input_sample
        dict
        fasta
        fasta_fai
        index_alignment
        intervals_and_num_intervals
        intervals_for_preprocessing
        known_sites_indels
        known_sites_indels_tbi

    main:

    // To gather all QC reports for MultiQC
    reports          = Channel.empty()
    versions         = Channel.empty()

    // PREPROCESSING

    if (params.step == 'mapping') {

        // STEP 0: QC & TRIM
        // Trim only with `--trim_fastq`
        // Additional options to be set up

        // UMI consensus calling
        if (params.umi_read_structure) {
            FASTQ_CREATE_UMI_CONSENSUS_FGBIO(
                input_fastq,
                fasta,
                fasta_fai,
                index_alignment,
                params.group_by_umi_strategy)

            bam_converted_from_fastq = FASTQ_CREATE_UMI_CONSENSUS_FGBIO.out.consensusbam.map{ meta, bam -> [ meta, bam, [] ] }

            // Convert back to fastq for further preprocessing
            // fasta are not needed when converting bam to fastq -> [ id:"fasta" ], []
            // No need for fasta.fai -> []
            interleave_input = false // Currently don't allow interleaved input
            CONVERT_FASTQ_UMI(
                bam_converted_from_fastq,
                [ [ id:"fasta" ], [] ], // fasta
                [ [ id:'null' ], [] ],  // fasta_fai
                interleave_input)

            reads_for_fastp = CONVERT_FASTQ_UMI.out.reads

            // Gather used softwares versions
            versions = versions.mix(CONVERT_FASTQ_UMI.out.versions)
            versions = versions.mix(FASTQ_CREATE_UMI_CONSENSUS_FGBIO.out.versions)
        } else {
            reads_for_fastp = input_fastq
        }

        // Trimming and/or splitting
        if (params.trim_fastq || params.split_fastq > 0) {

            save_trimmed_fail = false
            save_merged = false
            FASTP(
                reads_for_fastp,
                [], // we are not using any adapter fastas at the moment
                false, // we don't use discard_trimmed_pass at the moment
                save_trimmed_fail,
                save_merged
            )

            reports = reports.mix(FASTP.out.json.collect{ _meta, json -> json })
            reports = reports.mix(FASTP.out.html.collect{ _meta, html -> html })

            if (params.split_fastq) {
                reads_for_alignment = FASTP.out.reads.map{ meta, reads ->
                    def read_files = reads.sort(false) { a,b -> a.getName().tokenize('.')[0] <=> b.getName().tokenize('.')[0] }.collate(2)
                    [ meta + [ n_fastq: read_files.size() ], read_files ]
                }.transpose()
            } else reads_for_alignment = FASTP.out.reads

            versions = versions.mix(FASTP.out.versions)

        } else {
            reads_for_alignment = reads_for_fastp
        }

        // STEP 1: MAPPING READS TO REFERENCE GENOME
        // First, we must calculate number of lanes for each sample (meta.n_fastq)
        // This is needed to group reads from the same sample together using groupKey to avoid stalling the workflow
        // when reads from different samples are mixed together
        reads_for_alignment.map { meta, reads ->
                [ meta.subMap('patient', 'sample', 'sex', 'status'), reads ]
            }
            .groupTuple()
            .map { meta, reads ->
                meta + [ n_fastq: reads.size() ] // We can drop the FASTQ files now that we know how many there are
            }
            .set { reads_grouping_key }

        reads_for_alignment = reads_for_alignment.map{ meta, reads ->
            // Update meta.id to meta.sample no multiple lanes or splitted fastqs
            if (meta.size * meta.num_lanes == 1) [ meta + [ id:meta.sample ], reads ]
            else [ meta, reads ]
        }

        // reads will be sorted
        sort_bam = true
        FASTQ_ALIGN(reads_for_alignment, index_alignment, sort_bam, fasta, fasta_fai)

        // Grouping the bams from the same samples not to stall the workflow
        // Use groupKey to make sure that the correct group can advance as soon as it is complete
        // and not stall the workflow until all reads from all channels are mapped
        bam_mapped = FASTQ_ALIGN.out.bam
            .combine(reads_grouping_key) // Creates a tuple of [ meta, bam, reads_grouping_key ]
            .filter { meta1, _bam, meta2 -> meta1.sample == meta2.sample }
            // Add n_fastq and other variables to meta
            .map { meta1, bam, meta2 ->
                [ meta1 + meta2, bam ]
            }
            // Manipulate meta map to remove old fields and add new ones
            .map { meta, bam ->
                [ meta - meta.subMap('id', 'read_group', 'data_type', 'num_lanes', 'read_group', 'size') + [ data_type: 'bam', id: meta.sample ], bam ]
            }
            // Create groupKey from meta map
            .map { meta, bam ->
                [ groupKey( meta, meta.n_fastq), bam ]
            }
            // Group
            .groupTuple()

        bai_mapped = FASTQ_ALIGN.out.bai
            .combine(reads_grouping_key) // Creates a tuple of [ meta, bai, reads_grouping_key ]
            .filter { meta1, _bai, meta2 -> meta1.sample == meta2.sample }
            // Add n_fastq and other variables to meta
            .map { meta1, bai, meta2 ->
                [ meta1 + meta2, bai ]
            }
            // Manipulate meta map to remove old fields and add new ones
            .map { meta, bai ->
                [ meta - meta.subMap('id', 'read_group', 'data_type', 'num_lanes', 'read_group', 'size') + [ data_type: 'bai', id: meta.sample ], bai ]
            }
            // Create groupKey from meta map
            .map { meta, bai ->
                [ groupKey( meta, meta.n_fastq), bai ]
            }
            // Group
            .groupTuple()


        // gatk4 markduplicates can handle multiple bams as input, so no need to merge/index here
        // Except if and only if save_mapped or (skipping markduplicates and sentieon-dedup)
        if (
            params.save_mapped ||
            (
                (params.skip_tools && params.skip_tools.split(',').contains('markduplicates')) &&
                !(params.tools && params.tools.split(',').contains('sentieon_dedup'))
            )
        ) {
            // bams are merged (when multiple lanes from the same sample), indexed and then converted to cram
            BAM_MERGE_INDEX_SAMTOOLS(bam_mapped)

            BAM_TO_CRAM_MAPPING(BAM_MERGE_INDEX_SAMTOOLS.out.bam_bai, fasta, fasta_fai)
            // Create CSV to restart from this step
            if (params.save_output_as_bam) CHANNEL_ALIGN_CREATE_CSV(BAM_MERGE_INDEX_SAMTOOLS.out.bam_bai, params.outdir, params.save_output_as_bam)
            else CHANNEL_ALIGN_CREATE_CSV(BAM_TO_CRAM_MAPPING.out.cram.join(BAM_TO_CRAM_MAPPING.out.crai, failOnDuplicate: true, failOnMismatch: true), params.outdir, params.save_output_as_bam)

            // Gather used softwares versions
            versions = versions.mix(BAM_MERGE_INDEX_SAMTOOLS.out.versions)
            versions = versions.mix(BAM_TO_CRAM_MAPPING.out.versions)
        }

        // Gather used softwares versions
        versions = versions.mix(FASTQ_ALIGN.out.versions)
    }

    if (params.step in ['mapping', 'markduplicates']) {

        // ch_cram_no_markduplicates_restart = Channel.empty()
        cram_markduplicates_no_spark = Channel.empty()
        cram_sentieon_dedup          = Channel.empty()
        cram_markduplicates_spark    = Channel.empty()

        // STEP 2: markduplicates (+QC) + convert to CRAM

        // ch_bam_for_markduplicates will contain bam mapped with FASTQ_ALIGN when step is mapping
        // Or bams that are specified in the samplesheet.csv when step is prepare_recalibration
        cram_for_markduplicates = params.step == 'mapping' ? bam_mapped : input_sample.map{ meta, input, _index -> [ meta, input ] }
        // if no MD is done, then run QC on mapped & converted CRAM files
        // or the input BAM (+converted) or CRAM files
        cram_skip_markduplicates = Channel.empty()

        // Should it be possible to restart from converted crams?
        // For now, conversion from bam to cram is only done when skipping markduplicates

        if (
            params.skip_tools &&
            params.skip_tools.split(',').contains('markduplicates') &&
            !(params.tools && params.tools.split(',').contains('sentieon_dedup'))
        ) {
            if (params.step == 'mapping') {
                cram_skip_markduplicates = BAM_TO_CRAM_MAPPING.out.cram.join(BAM_TO_CRAM_MAPPING.out.crai, failOnDuplicate: true, failOnMismatch: true)
            } else {
                cram_skip_markduplicates = Channel.empty().mix(input_sample)
            }

            CRAM_QC_NO_MD(cram_skip_markduplicates, fasta, intervals_for_preprocessing)

            // Gather QC reports
            reports = reports.mix(CRAM_QC_NO_MD.out.reports.collect{ _meta, report -> [ report ] })

            // Gather used softwares versions
            versions = versions.mix(CRAM_QC_NO_MD.out.versions)
        } else if (params.use_gatk_spark && params.use_gatk_spark.contains('markduplicates')) {
            BAM_MARKDUPLICATES_SPARK(
                cram_for_markduplicates,
                dict,
                fasta,
                fasta_fai,
                intervals_for_preprocessing)
            cram_markduplicates_spark = BAM_MARKDUPLICATES_SPARK.out.cram

            // Gather QC reports
            reports = reports.mix(BAM_MARKDUPLICATES_SPARK.out.reports.collect{ _meta, report -> [ report ] })

            // Gather used softwares versions
            versions = versions.mix(BAM_MARKDUPLICATES_SPARK.out.versions)
        } else if (params.tools && params.tools.split(',').contains('sentieon_dedup')) {
            crai_for_markduplicates = params.step == 'mapping' ? bai_mapped : input_sample.map{ meta, _input, index -> [ meta, index ] }
            BAM_SENTIEON_DEDUP(
                cram_for_markduplicates,
                crai_for_markduplicates,
                fasta,
                fasta_fai,
                intervals_for_preprocessing)

            cram_sentieon_dedup = BAM_SENTIEON_DEDUP.out.cram

            // Gather QC reports
            reports = reports.mix(BAM_SENTIEON_DEDUP.out.reports.collect{ _meta, report -> [ report ] })

            // Gather used softwares versions
            versions = versions.mix(BAM_SENTIEON_DEDUP.out.versions)
        } else {
            BAM_MARKDUPLICATES(
                cram_for_markduplicates,
                fasta,
                fasta_fai,
                intervals_for_preprocessing)

            cram_markduplicates_no_spark = BAM_MARKDUPLICATES.out.cram

            // Gather QC reports
            reports = reports.mix(BAM_MARKDUPLICATES.out.reports.collect{ _meta, report -> [ report ] })

            // Gather used softwares versions
            versions = versions.mix(BAM_MARKDUPLICATES.out.versions)
        }

        // ch_md_cram_for_restart contains either:
        // - crams from markduplicates
        // - crams from sentieon_dedup
        // - crams from markduplicates_spark
        // - crams from input step markduplicates --> from the converted ones only?
        ch_md_cram_for_restart = Channel.empty().mix(cram_markduplicates_no_spark, cram_markduplicates_spark, cram_sentieon_dedup)
            // Make sure correct data types are carried through
            .map{ meta, cram, crai -> [ meta + [data_type: "cram"], cram, crai ] }

        // If params.save_output_as_bam, then convert CRAM files to BAM
        CRAM_TO_BAM(ch_md_cram_for_restart, fasta, fasta_fai)
        versions = versions.mix(CRAM_TO_BAM.out.versions)

        // CSV should be written for the file actually out, either CRAM or BAM
        // Create CSV to restart from this step
        csv_subfolder = (params.tools && params.tools.split(',').contains('sentieon_dedup')) ? 'sentieon_dedup' : 'markduplicates'

        if (params.save_output_as_bam) CHANNEL_MARKDUPLICATES_CREATE_CSV(CRAM_TO_BAM.out.bam.join(CRAM_TO_BAM.out.bai, failOnDuplicate: true, failOnMismatch: true), csv_subfolder, params.outdir, params.save_output_as_bam)
        else CHANNEL_MARKDUPLICATES_CREATE_CSV(ch_md_cram_for_restart, csv_subfolder, params.outdir, params.save_output_as_bam)
    }

    if (params.step in ['mapping', 'markduplicates', 'prepare_recalibration']) {

        // Run if starting from step "prepare_recalibration"
        if (params.step == 'prepare_recalibration') {

            ch_cram_for_bam_baserecalibrator = Channel.empty().mix(input_sample)

            // Set the input samples for restart so we generate a samplesheet that contains the input files together with the recalibration table
            ch_md_cram_for_restart           = ch_cram_for_bam_baserecalibrator

        } else {

            // ch_cram_for_bam_baserecalibrator contains either:
            // - crams from markduplicates
            // - crams from markduplicates_spark
            // - crams converted from bam mapped when skipping markduplicates
            // - input cram files, when start from step markduplicates
            ch_cram_for_bam_baserecalibrator = Channel.empty().mix(ch_md_cram_for_restart, cram_skip_markduplicates )
                // Make sure correct data types are carried through
                .map{ meta, cram, crai -> [ meta + [data_type: "cram"], cram, crai ] }

        }

        // STEP 3: Create recalibration tables
        if (!(params.skip_tools && params.skip_tools.split(',').contains('baserecalibrator'))) {

            ch_table_bqsr_no_spark = Channel.empty()
            ch_table_bqsr_spark    = Channel.empty()

            if (params.use_gatk_spark && params.use_gatk_spark.contains('baserecalibrator')) {
            BAM_BASERECALIBRATOR_SPARK(
                ch_cram_for_bam_baserecalibrator,
                dict,
                fasta,
                fasta_fai,
                intervals_and_num_intervals,
                known_sites_indels,
                known_sites_indels_tbi)

                ch_table_bqsr_spark = BAM_BASERECALIBRATOR_SPARK.out.table_bqsr

                // Gather used softwares versions
                versions = versions.mix(BAM_BASERECALIBRATOR_SPARK.out.versions)
            } else {

            BAM_BASERECALIBRATOR(
                ch_cram_for_bam_baserecalibrator,
                dict,
                fasta,
                fasta_fai,
                intervals_and_num_intervals,
                known_sites_indels,
                known_sites_indels_tbi)

                ch_table_bqsr_no_spark = BAM_BASERECALIBRATOR.out.table_bqsr

                // Gather used softwares versions
                versions = versions.mix(BAM_BASERECALIBRATOR.out.versions)
            }

            // ch_table_bqsr contains either:
            // - bqsr table from baserecalibrator
            // - bqsr table from baserecalibrator_spark
            ch_table_bqsr = Channel.empty().mix(
                ch_table_bqsr_no_spark,
                ch_table_bqsr_spark)

            reports = reports.mix(ch_table_bqsr.collect{ _meta, table -> [ table ] })

            cram_applybqsr = ch_cram_for_bam_baserecalibrator.join(ch_table_bqsr, failOnDuplicate: true, failOnMismatch: true)

            // Create CSV to restart from this step
            CHANNEL_BASERECALIBRATOR_CREATE_CSV(ch_md_cram_for_restart.join(ch_table_bqsr, failOnDuplicate: true), params.tools, params.skip_tools, params.outdir, params.save_output_as_bam)
        }
    }

    // STEP 4: RECALIBRATING
    if (params.step in ['mapping', 'markduplicates', 'prepare_recalibration', 'recalibrate']) {

        // Run if starting from step "prepare_recalibration"
        if (params.step == 'recalibrate') {

            cram_applybqsr = Channel.empty().mix(input_sample)

        }
        cram_applybqsr.view()
        if (!(params.skip_tools && params.skip_tools.split(',').contains('baserecalibrator'))) {
            cram_variant_calling_no_spark = Channel.empty()
            cram_variant_calling_spark    = Channel.empty()

            if (params.use_gatk_spark && params.use_gatk_spark.contains('baserecalibrator')) {

                BAM_APPLYBQSR_SPARK(
                    cram_applybqsr,
                    dict,
                    fasta,
                    fasta_fai,
                    intervals_and_num_intervals)

                cram_variant_calling_spark = BAM_APPLYBQSR_SPARK.out.cram

                // Gather used softwares versions
                versions = versions.mix(BAM_APPLYBQSR_SPARK.out.versions)

            } else {

                BAM_APPLYBQSR(
                    cram_applybqsr,
                    dict,
                    fasta,
                    fasta_fai,
                    intervals_and_num_intervals)

                cram_variant_calling_no_spark = BAM_APPLYBQSR.out.cram

                // Gather used softwares versions
                versions = versions.mix(BAM_APPLYBQSR.out.versions)
            }

            cram_variant_calling = Channel.empty().mix(
                cram_variant_calling_no_spark,
                cram_variant_calling_spark)

            // If params.save_output_as_bam, then convert CRAM files to BAM
            CRAM_TO_BAM_RECAL(cram_variant_calling, fasta, fasta_fai)
            versions = versions.mix(CRAM_TO_BAM_RECAL.out.versions)

            // CSV should be written for the file actually out out, either CRAM or BAM
            csv_recalibration = Channel.empty()
            csv_recalibration = params.save_output_as_bam ? CRAM_TO_BAM_RECAL.out.bam.join(CRAM_TO_BAM_RECAL.out.bai, failOnDuplicate: true, failOnMismatch: true) : cram_variant_calling

            // Create CSV to restart from this step
            CHANNEL_APPLYBQSR_CREATE_CSV(csv_recalibration, params.outdir, params.save_output_as_bam)

        } else if (params.step == 'recalibrate') {
            // cram_variant_calling contains either:
            // - input bams converted to crams, if started from step recal + skip BQSR
            // - input crams if started from step recal + skip BQSR
            cram_variant_calling = Channel.empty().mix(input_sample.map{ meta, cram, crai, _table -> [ meta, cram, crai ] })
        } else {
            // cram_variant_calling contains either:
            // - crams from markduplicates = ch_cram_for_bam_baserecalibrator if skip BQSR but not started from step recalibration
            cram_variant_calling = Channel.empty().mix(ch_cram_for_bam_baserecalibrator)
        }
    }

    emit:
    cram_variant_calling
    reports
    versions

}
