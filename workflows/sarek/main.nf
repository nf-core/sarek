/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { paramsSummaryMap                                  } from 'plugin/nf-schema'
include { paramsSummaryMultiqc                              } from '../../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML                            } from '../../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText                            } from '../../subworkflows/local/utils_nfcore_sarek_pipeline'

// Create samplesheets to restart from different steps
include { CHANNEL_ALIGN_CREATE_CSV                          } from '../../subworkflows/local/channel_align_create_csv/main'
include { CHANNEL_MARKDUPLICATES_CREATE_CSV                 } from '../../subworkflows/local/channel_markduplicates_create_csv/main'
include { CHANNEL_BASERECALIBRATOR_CREATE_CSV               } from '../../subworkflows/local/channel_baserecalibrator_create_csv/main'
include { CHANNEL_APPLYBQSR_CREATE_CSV                      } from '../../subworkflows/local/channel_applybqsr_create_csv/main'
include { CHANNEL_VARIANT_CALLING_CREATE_CSV                } from '../../subworkflows/local/channel_variant_calling_create_csv/main'


// Convert BAM files to FASTQ files
include { BAM_CONVERT_SAMTOOLS as CONVERT_FASTQ_INPUT       } from '../../subworkflows/local/bam_convert_samtools/main'
include { BAM_CONVERT_SAMTOOLS as CONVERT_FASTQ_UMI         } from '../../subworkflows/local/bam_convert_samtools/main'

// Convert fastq.gz.spring files to fastq.gz files
include { SPRING_DECOMPRESS as SPRING_DECOMPRESS_TO_R1_FQ   } from '../../modules/nf-core/spring/decompress/main'
include { SPRING_DECOMPRESS as SPRING_DECOMPRESS_TO_R2_FQ   } from '../../modules/nf-core/spring/decompress/main'
include { SPRING_DECOMPRESS as SPRING_DECOMPRESS_TO_FQ_PAIR } from '../../modules/nf-core/spring/decompress/main'

// Run FASTQC
include { FASTQC                                            } from '../../modules/nf-core/fastqc/main'

// TRIM/SPLIT FASTQ Files
include { FASTP                                             } from '../../modules/nf-core/fastp/main'

// Create umi consensus bams from fastq
include { FASTQ_CREATE_UMI_CONSENSUS_FGBIO                  } from '../../subworkflows/local/fastq_create_umi_consensus_fgbio/main'

// Map input reads to reference genome
include { FASTQ_ALIGN_BWAMEM_MEM2_DRAGMAP_SENTIEON          } from '../../subworkflows/local/fastq_align_bwamem_mem2_dragmap_sentieon/main'

// Merge and index BAM files (optional)
include { BAM_MERGE_INDEX_SAMTOOLS                          } from '../../subworkflows/local/bam_merge_index_samtools/main'

// Convert BAM files
include { SAMTOOLS_CONVERT as BAM_TO_CRAM_MAPPING           } from '../../modules/nf-core/samtools/convert/main'

// Convert CRAM files (optional)
include { SAMTOOLS_CONVERT as CRAM_TO_BAM                   } from '../../modules/nf-core/samtools/convert/main'
include { SAMTOOLS_CONVERT as CRAM_TO_BAM_RECAL             } from '../../modules/nf-core/samtools/convert/main'

// Mark Duplicates (+QC)
include { BAM_MARKDUPLICATES                                } from '../../subworkflows/local/bam_markduplicates/main'
include { BAM_MARKDUPLICATES_SPARK                          } from '../../subworkflows/local/bam_markduplicates_spark/main'
include { BAM_SENTIEON_DEDUP                                } from '../../subworkflows/local/bam_sentieon_dedup/main'

// QC on CRAM
include { CRAM_QC_MOSDEPTH_SAMTOOLS as CRAM_QC_NO_MD        } from '../../subworkflows/local/cram_qc_mosdepth_samtools/main'
include { CRAM_SAMPLEQC                                     } from '../../subworkflows/local/cram_sampleqc/main'

// Create recalibration tables
include { BAM_BASERECALIBRATOR                              } from '../../subworkflows/local/bam_baserecalibrator/main'
include { BAM_BASERECALIBRATOR_SPARK                        } from '../../subworkflows/local/bam_baserecalibrator_spark/main'

// Create recalibrated cram files to use for variant calling (+QC)
include { BAM_APPLYBQSR                                     } from '../../subworkflows/local/bam_applybqsr/main'
include { BAM_APPLYBQSR_SPARK                               } from '../../subworkflows/local/bam_applybqsr_spark/main'

// Variant calling on a single normal sample
include { BAM_VARIANT_CALLING_GERMLINE_ALL                  } from '../../subworkflows/local/bam_variant_calling_germline_all/main'

// Variant calling on a single tumor sample
include { BAM_VARIANT_CALLING_TUMOR_ONLY_ALL                } from '../../subworkflows/local/bam_variant_calling_tumor_only_all/main'

// Variant calling on tumor/normal pair
include { BAM_VARIANT_CALLING_SOMATIC_ALL                   } from '../../subworkflows/local/bam_variant_calling_somatic_all/main'

// POST VARIANTCALLING: e.g. merging
include { POST_VARIANTCALLING                               } from '../../subworkflows/local/post_variantcalling/main'

// QC on VCF files
include { VCF_QC_BCFTOOLS_VCFTOOLS                          } from '../../subworkflows/local/vcf_qc_bcftools_vcftools/main'

// Annotation
include { VCF_ANNOTATE_ALL                                  } from '../../subworkflows/local/vcf_annotate_all/main'

// MULTIQC
include { MULTIQC                                           } from '../../modules/nf-core/multiqc/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow SAREK {
    take:
        input_sample
        allele_files
        bcftools_annotations
        bcftools_annotations_tbi
        bcftools_header_lines
        cf_chrom_len
        chr_files
        cnvkit_reference
        dbsnp
        dbsnp_tbi
        dbsnp_vqsr
        dict
        fasta
        fasta_fai
        gc_file
        germline_resource
        germline_resource_tbi
        index_alignment
        intervals_and_num_intervals
        intervals_bed_combined
        intervals_bed_combined_for_variant_calling
        intervals_bed_gz_tbi_and_num_intervals
        intervals_bed_gz_tbi_combined
        intervals_for_preprocessing
        known_indels_vqsr
        known_sites_indels
        known_sites_indels_tbi
        known_sites_snps
        known_sites_snps_tbi
        known_snps_vqsr
        loci_files
        mappability
        msisensorpro_scan
        ngscheckmate_bed
        pon
        pon_tbi
        rt_file
        sentieon_dnascope_model
        snpeff_cache
        vep_cache
        vep_cache_version
        vep_extra_files
        vep_fasta
        vep_genome
        vep_species

    main:

    // To gather all QC reports for MultiQC
    ch_multiqc_files = Channel.empty()
    multiqc_report   = Channel.empty()
    reports          = Channel.empty()
    versions         = Channel.empty()

    // PREPROCESSING

    if (params.step == 'mapping') {

        // Figure out if input is bam, fastq, or spring
        input_sample_type = input_sample.branch{
            bam:                 it[0].data_type == "bam"
            fastq_gz:            it[0].data_type == "fastq_gz"
            one_fastq_gz_spring: it[0].data_type == "one_fastq_gz_spring"
            two_fastq_gz_spring: it[0].data_type == "two_fastq_gz_spring"
        }

        // Two fastq.gz-files
        fastq_gz = input_sample_type.fastq_gz.map { meta, files -> addReadgroupToMeta(meta, files) }

        // Just one fastq.gz.spring-file with both R1 and R2
        fastq_gz_pair_from_spring = SPRING_DECOMPRESS_TO_FQ_PAIR(input_sample_type.one_fastq_gz_spring, false)

        one_fastq_gz_from_spring = fastq_gz_pair_from_spring.fastq.map { meta, files -> addReadgroupToMeta(meta, files) }

        // Two fastq.gz.spring-files - one for R1 and one for R2
        r1_fastq_gz_from_spring = SPRING_DECOMPRESS_TO_R1_FQ(input_sample_type.two_fastq_gz_spring.map{ meta, files ->
            [meta, files[0] ]},
            true // write_one_fastq_gz
        )
        r2_fastq_gz_from_spring = SPRING_DECOMPRESS_TO_R2_FQ(input_sample_type.two_fastq_gz_spring.map{ meta, files ->
            [meta, files[1] ]},
            true // write_one_fastq_gz
        )

        two_fastq_gz_from_spring = r1_fastq_gz_from_spring.fastq.join(r2_fastq_gz_from_spring.fastq).map{ meta, fastq_1, fastq_2 -> [meta, [fastq_1, fastq_2]]}

        two_fastq_gz_from_spring = two_fastq_gz_from_spring.map { meta, files -> addReadgroupToMeta(meta, files) }

        // Convert any bam input to fastq
        // fasta are not needed when converting bam to fastq -> [ id:"fasta" ], []
        // No need for fasta.fai -> []
        interleave_input = false // Currently don't allow interleaved input
        CONVERT_FASTQ_INPUT(
            input_sample_type.bam,
            [ [ id:"fasta" ], [] ], // fasta
            [ [ id:'null' ], [] ],  // fasta_fai
            interleave_input)

        // Gather fastq (inputed or converted)
        // Theorically this could work on mixed input (fastq for one sample and bam for another)
        // But not sure how to handle that with the samplesheet
        // Or if we really want users to be able to do that
        input_fastq = fastq_gz.mix(CONVERT_FASTQ_INPUT.out.reads).mix(one_fastq_gz_from_spring).mix(two_fastq_gz_from_spring)

        // STEP 0: QC & TRIM
        // `--skip_tools fastqc` to skip fastqc
        // Trim only with `--trim_fastq`
        // Additional options to be set up

        // QC
        if (!(params.skip_tools && params.skip_tools.split(',').contains('fastqc'))) {
            FASTQC(input_fastq)

            reports = reports.mix(FASTQC.out.zip.collect{ meta, logs -> logs })
            versions = versions.mix(FASTQC.out.versions.first())
        }

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

            reports = reports.mix(FASTP.out.json.collect{ meta, json -> json })
            reports = reports.mix(FASTP.out.html.collect{ meta, html -> html })

            if (params.split_fastq) {
                reads_for_alignment = FASTP.out.reads.map{ meta, reads ->
                    read_files = reads.sort(false) { a,b -> a.getName().tokenize('.')[0] <=> b.getName().tokenize('.')[0] }.collate(2)
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
        FASTQ_ALIGN_BWAMEM_MEM2_DRAGMAP_SENTIEON(reads_for_alignment, index_alignment, sort_bam, fasta, fasta_fai)

        // Grouping the bams from the same samples not to stall the workflow
        // Use groupKey to make sure that the correct group can advance as soon as it is complete
        // and not stall the workflow until all reads from all channels are mapped
        bam_mapped = FASTQ_ALIGN_BWAMEM_MEM2_DRAGMAP_SENTIEON.out.bam
            .combine(reads_grouping_key) // Creates a tuple of [ meta, bam, reads_grouping_key ]
            .filter { meta1, bam, meta2 -> meta1.sample == meta2.sample }
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

        bai_mapped = FASTQ_ALIGN_BWAMEM_MEM2_DRAGMAP_SENTIEON.out.bai
            .combine(reads_grouping_key) // Creates a tuple of [ meta, bai, reads_grouping_key ]
            .filter { meta1, bai, meta2 -> meta1.sample == meta2.sample }
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
        versions = versions.mix(CONVERT_FASTQ_INPUT.out.versions)
        versions = versions.mix(FASTQ_ALIGN_BWAMEM_MEM2_DRAGMAP_SENTIEON.out.versions)
    }

    if (params.step in ['mapping', 'markduplicates']) {

        // ch_cram_no_markduplicates_restart = Channel.empty()
        cram_markduplicates_no_spark = Channel.empty()
        cram_sentieon_dedup          = Channel.empty()
        cram_markduplicates_spark    = Channel.empty()

        // STEP 2: markduplicates (+QC) + convert to CRAM

        // ch_bam_for_markduplicates will contain bam mapped with FASTQ_ALIGN_BWAMEM_MEM2_DRAGMAP_SENTIEON when step is mapping
        // Or bams that are specified in the samplesheet.csv when step is prepare_recalibration
        cram_for_markduplicates = params.step == 'mapping' ? bam_mapped : input_sample.map{ meta, input, index -> [ meta, input ] }
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
            reports = reports.mix(CRAM_QC_NO_MD.out.reports.collect{ meta, report -> [ report ] })

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
            reports = reports.mix(BAM_MARKDUPLICATES_SPARK.out.reports.collect{ meta, report -> [ report ] })

            // Gather used softwares versions
            versions = versions.mix(BAM_MARKDUPLICATES_SPARK.out.versions)
        } else if (params.tools && params.tools.split(',').contains('sentieon_dedup')) {
            crai_for_markduplicates = params.step == 'mapping' ? bai_mapped : input_sample.map{ meta, input, index -> [ meta, index ] }
            BAM_SENTIEON_DEDUP(
                cram_for_markduplicates,
                crai_for_markduplicates,
                fasta,
                fasta_fai,
                intervals_for_preprocessing)

            cram_sentieon_dedup = BAM_SENTIEON_DEDUP.out.cram

            // Gather QC reports
            reports = reports.mix(BAM_SENTIEON_DEDUP.out.reports.collect{ meta, report -> [ report ] })

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
            reports = reports.mix(BAM_MARKDUPLICATES.out.reports.collect{ meta, report -> [ report ] })

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

            reports = reports.mix(ch_table_bqsr.collect{ meta, table -> [ table ] })

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
            cram_variant_calling = Channel.empty().mix(input_sample.map{ meta, cram, crai, table -> [ meta, cram, crai ] })
        } else {
            // cram_variant_calling contains either:
            // - crams from markduplicates = ch_cram_for_bam_baserecalibrator if skip BQSR but not started from step recalibration
            cram_variant_calling = Channel.empty().mix(ch_cram_for_bam_baserecalibrator)
        }
    }


    if (params.step == 'variant_calling') {

        cram_variant_calling = Channel.empty().mix( input_sample )

    }

    if (params.step == 'annotate') cram_variant_calling = Channel.empty()

        // RUN CRAM QC on the recalibrated CRAM files or when starting from step variant calling. NGSCheckmate should be run also on non-recalibrated CRAM files
        CRAM_SAMPLEQC(cram_variant_calling,
            ngscheckmate_bed,
            fasta,
            params.skip_tools && params.skip_tools.split(',').contains('baserecalibrator'),
            intervals_for_preprocessing)

    if (params.tools) {

        //
        // Logic to separate germline samples, tumor samples with no matched normal, and combine tumor-normal pairs
        //
        cram_variant_calling_status = cram_variant_calling.branch{
            normal: it[0].status == 0
            tumor:  it[0].status == 1
        }

        // All Germline samples
        cram_variant_calling_normal_to_cross = cram_variant_calling_status.normal.map{ meta, cram, crai -> [ meta.patient, meta, cram, crai ] }

        // All tumor samples
        cram_variant_calling_pair_to_cross = cram_variant_calling_status.tumor.map{ meta, cram, crai -> [ meta.patient, meta, cram, crai ] }

        // Tumor only samples
        // 1. Group together all tumor samples by patient ID [ patient1, [ meta1, meta2 ], [ cram1, crai1, cram2, crai2 ] ]

        // Downside: this only works by waiting for all tumor samples to finish preprocessing, since no group size is provided
        cram_variant_calling_tumor_grouped = cram_variant_calling_pair_to_cross.groupTuple()

        // 2. Join with normal samples, in each channel there is one key per patient now. Patients without matched normal end up with: [ patient1, [ meta1, meta2 ], [ cram1, crai1, cram2, crai2 ], null ]
        cram_variant_calling_tumor_joined = cram_variant_calling_tumor_grouped.join(cram_variant_calling_normal_to_cross, failOnDuplicate: true, remainder: true)

        // 3. Filter out entries with last entry null
        cram_variant_calling_tumor_filtered = cram_variant_calling_tumor_joined.filter{ it ->  !(it.last()) }

        // 4. Transpose [ patient1, [ meta1, meta2 ], [ cram1, crai1, cram2, crai2 ] ] back to [ patient1, meta1, [ cram1, crai1 ], null ] [ patient1, meta2, [ cram2, crai2 ], null ]
        // and remove patient ID field & null value for further processing [ meta1, [ cram1, crai1 ] ] [ meta2, [ cram2, crai2 ] ]
        cram_variant_calling_tumor_only = cram_variant_calling_tumor_filtered.transpose().map{ it -> [it[1], it[2], it[3]] }

        if (params.only_paired_variant_calling) {
            // Normal only samples

            // 1. Join with tumor samples, in each channel there is one key per patient now. Patients without matched tumor end up with: [ patient1, [ meta1 ], [ cram1, crai1 ], null ] as there is only one matched normal possible
            cram_variant_calling_normal_joined = cram_variant_calling_normal_to_cross.join(cram_variant_calling_tumor_grouped, failOnDuplicate: true, remainder: true)

            // 2. Filter out entries with last entry null
            cram_variant_calling_normal_filtered = cram_variant_calling_normal_joined.filter{ it ->  !(it.last()) }

            // 3. Remove patient ID field & null value for further processing [ meta1, [ cram1, crai1 ] ] [ meta2, [ cram2, crai2 ] ] (no transposing needed since only one normal per patient ID)
            cram_variant_calling_status_normal = cram_variant_calling_normal_filtered.map{ it -> [it[1], it[2], it[3]] }

        } else {
            cram_variant_calling_status_normal = cram_variant_calling_status.normal
        }

        // Tumor - normal pairs
        // Use cross to combine normal with all tumor samples, i.e. multi tumor samples from recurrences
        cram_variant_calling_pair = cram_variant_calling_normal_to_cross.cross(cram_variant_calling_pair_to_cross)
            .map { normal, tumor ->
                def meta = [:]

                meta.id         = "${tumor[1].sample}_vs_${normal[1].sample}".toString()
                meta.normal_id  = normal[1].sample
                meta.patient    = normal[0]
                meta.sex        = normal[1].sex
                meta.tumor_id   = tumor[1].sample

                [ meta, normal[2], normal[3], tumor[2], tumor[3] ]
            }

        // GERMLINE VARIANT CALLING
        BAM_VARIANT_CALLING_GERMLINE_ALL(
            params.tools,
            params.skip_tools,
            cram_variant_calling_status_normal,
            [ [ id:'bwa' ], [] ], // bwa_index for tiddit; not used here
            cnvkit_reference,
            dbsnp,
            dbsnp_tbi,
            dbsnp_vqsr,
            dict,
            fasta,
            fasta_fai,
            intervals_and_num_intervals,
            intervals_bed_combined, // [] if no_intervals, else interval_bed_combined.bed,
            intervals_bed_gz_tbi_combined, // [] if no_intervals, else interval_bed_combined_gz, interval_bed_combined_gz_tbi
            intervals_bed_combined_for_variant_calling, // no_intervals.bed if no intervals, else interval_bed_combined.bed; Channel operations possible
            intervals_bed_gz_tbi_and_num_intervals,
            known_indels_vqsr,
            known_sites_indels,
            known_sites_indels_tbi,
            known_sites_snps,
            known_sites_snps_tbi,
            known_snps_vqsr,
            params.joint_germline,
            params.skip_tools && params.skip_tools.split(',').contains('haplotypecaller_filter'), // true if filtering should be skipped
            params.sentieon_haplotyper_emit_mode,
            params.sentieon_dnascope_emit_mode,
            params.sentieon_dnascope_pcr_indel_model,
            sentieon_dnascope_model)

        // TUMOR ONLY VARIANT CALLING
        BAM_VARIANT_CALLING_TUMOR_ONLY_ALL(
            params.tools,
            cram_variant_calling_tumor_only,
            [ [ id:'bwa' ], [] ], // bwa_index for tiddit; not used here
            cf_chrom_len,
            chr_files,
            cnvkit_reference,
            dbsnp,
            dbsnp_tbi,
            dict,
            fasta,
            fasta_fai,
            germline_resource,
            germline_resource_tbi,
            intervals_and_num_intervals,
            intervals_bed_gz_tbi_and_num_intervals,
            intervals_bed_combined,
            intervals_bed_gz_tbi_combined, // [] if no_intervals, else interval_bed_combined_gz, interval_bed_combined_gz_tbi
            mappability,
            pon,
            pon_tbi,
            params.joint_mutect2,
            params.wes
        )

        // PAIR VARIANT CALLING
        BAM_VARIANT_CALLING_SOMATIC_ALL(
            params.tools,
            cram_variant_calling_pair,
            [ [ id:'bwa' ], [] ], // bwa_index for tiddit; not used here
            cf_chrom_len,
            chr_files,
            dbsnp,
            dbsnp_tbi,
            dict,
            fasta,
            fasta_fai,
            germline_resource,
            germline_resource_tbi,
            intervals_and_num_intervals,
            intervals_bed_gz_tbi_and_num_intervals,
            intervals_bed_combined,
            intervals_bed_gz_tbi_combined, // [] if no_intervals, else interval_bed_combined_gz, interval_bed_combined_gz_tbi
            mappability,
            msisensorpro_scan,
            pon,
            pon_tbi,
            allele_files,
            loci_files,
            gc_file,
            rt_file,
            params.joint_mutect2,
            params.wes
        )

        // POST VARIANTCALLING
        POST_VARIANTCALLING(BAM_VARIANT_CALLING_GERMLINE_ALL.out.vcf_all,
                            params.concatenate_vcfs)

        // Gather vcf files for annotation and QC
        vcf_to_annotate = Channel.empty()
        vcf_to_annotate = vcf_to_annotate.mix(BAM_VARIANT_CALLING_GERMLINE_ALL.out.vcf_deepvariant)
        vcf_to_annotate = vcf_to_annotate.mix(BAM_VARIANT_CALLING_GERMLINE_ALL.out.vcf_freebayes)
        vcf_to_annotate = vcf_to_annotate.mix(BAM_VARIANT_CALLING_GERMLINE_ALL.out.vcf_haplotypecaller)
        vcf_to_annotate = vcf_to_annotate.mix(BAM_VARIANT_CALLING_GERMLINE_ALL.out.vcf_manta)
        vcf_to_annotate = vcf_to_annotate.mix(BAM_VARIANT_CALLING_GERMLINE_ALL.out.vcf_sentieon_dnascope)
        vcf_to_annotate = vcf_to_annotate.mix(BAM_VARIANT_CALLING_GERMLINE_ALL.out.vcf_sentieon_haplotyper)
        vcf_to_annotate = vcf_to_annotate.mix(BAM_VARIANT_CALLING_GERMLINE_ALL.out.vcf_strelka)
        vcf_to_annotate = vcf_to_annotate.mix(BAM_VARIANT_CALLING_GERMLINE_ALL.out.vcf_tiddit)
        vcf_to_annotate = vcf_to_annotate.mix(BAM_VARIANT_CALLING_GERMLINE_ALL.out.vcf_mpileup)
        vcf_to_annotate = vcf_to_annotate.mix(BAM_VARIANT_CALLING_TUMOR_ONLY_ALL.out.vcf_all)
        vcf_to_annotate = vcf_to_annotate.mix(BAM_VARIANT_CALLING_SOMATIC_ALL.out.vcf_all)

        // QC
        VCF_QC_BCFTOOLS_VCFTOOLS(vcf_to_annotate, intervals_bed_combined)

        reports = reports.mix(VCF_QC_BCFTOOLS_VCFTOOLS.out.bcftools_stats.collect{ meta, stats -> [ stats ] })
        reports = reports.mix(VCF_QC_BCFTOOLS_VCFTOOLS.out.vcftools_tstv_counts.collect{ meta, counts -> [ counts ] })
        reports = reports.mix(VCF_QC_BCFTOOLS_VCFTOOLS.out.vcftools_tstv_qual.collect{ meta, qual -> [ qual ] })
        reports = reports.mix(VCF_QC_BCFTOOLS_VCFTOOLS.out.vcftools_filter_summary.collect{ meta, summary -> [ summary ] })
        reports = reports.mix(BAM_VARIANT_CALLING_GERMLINE_ALL.out.out_indexcov.collect{ meta, indexcov -> indexcov.flatten() })
        reports = reports.mix(BAM_VARIANT_CALLING_SOMATIC_ALL.out.out_indexcov.collect{ meta, indexcov -> indexcov.flatten() })

        CHANNEL_VARIANT_CALLING_CREATE_CSV(vcf_to_annotate, params.outdir)

        // Gather used variant calling softwares versions
        versions = versions.mix(BAM_VARIANT_CALLING_GERMLINE_ALL.out.versions)
        versions = versions.mix(BAM_VARIANT_CALLING_SOMATIC_ALL.out.versions)
        versions = versions.mix(BAM_VARIANT_CALLING_TUMOR_ONLY_ALL.out.versions)
        versions = versions.mix(POST_VARIANTCALLING.out.versions)
        versions = versions.mix(VCF_QC_BCFTOOLS_VCFTOOLS.out.versions)

        // ANNOTATE
        if (params.step == 'annotate') vcf_to_annotate = input_sample

        if (params.tools.split(',').contains('merge') || params.tools.split(',').contains('snpeff') || params.tools.split(',').contains('vep')|| params.tools.split(',').contains('bcfann')) {

            vep_fasta = (params.vep_include_fasta) ? fasta : [[id: 'null'], []]

            VCF_ANNOTATE_ALL(
                vcf_to_annotate.map{meta, vcf -> [ meta + [ file_name: vcf.baseName ], vcf ] },
                vep_fasta,
                params.tools,
                params.snpeff_db,
                snpeff_cache,
                vep_genome,
                vep_species,
                vep_cache_version,
                vep_cache,
                vep_extra_files,
                bcftools_annotations,
                bcftools_annotations_tbi,
                bcftools_header_lines)

            // Gather used softwares versions
            versions = versions.mix(VCF_ANNOTATE_ALL.out.versions)
            reports = reports.mix(VCF_ANNOTATE_ALL.out.reports)
        }
    }

    //
    // Collate and save software versions
    //
    version_yaml = Channel.empty()
    if (!(params.skip_tools && params.skip_tools.split(',').contains('versions'))) {
        version_yaml = softwareVersionsToYAML(versions)
            .collectFile(storeDir: "${params.outdir}/pipeline_info", name: 'nf_core_sarek_software_mqc_versions.yml', sort: true, newLine: true)
    }

    //
    // MODULE: MultiQC
    //
    if (!(params.skip_tools && params.skip_tools.split(',').contains('multiqc'))) {

        ch_multiqc_config                     = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
        ch_multiqc_custom_config              = params.multiqc_config ? Channel.fromPath(params.multiqc_config, checkIfExists: true) : Channel.empty()
        ch_multiqc_logo                       = params.multiqc_logo ? Channel.fromPath(params.multiqc_logo, checkIfExists: true) : Channel.empty()
        summary_params                        = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
        ch_workflow_summary                   = Channel.value(paramsSummaryMultiqc(summary_params))
        ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
        ch_methods_description                = Channel.value(methodsDescriptionText(ch_multiqc_custom_methods_description))
        ch_multiqc_files                      = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
        ch_multiqc_files                      = ch_multiqc_files.mix(version_yaml)
        ch_multiqc_files                      = ch_multiqc_files.mix(reports)
        ch_multiqc_files                      = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml', sort: true))

        MULTIQC (
            ch_multiqc_files.collect(),
            ch_multiqc_config.toList(),
            ch_multiqc_custom_config.toList(),
            ch_multiqc_logo.toList(),
            [],
            []
        )
        multiqc_report = MULTIQC.out.report.toList()

    }

    emit:
    multiqc_report // channel: /path/to/multiqc_report.html
    versions       // channel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Add readgroup to meta and remove lane
def addReadgroupToMeta(meta, files) {
    def CN = params.seq_center ? "CN:${params.seq_center}\\t" : ''
    def flowcell = flowcellLaneFromFastq(files[0])

    // Check if flowcell ID matches
    if ( flowcell && flowcell != flowcellLaneFromFastq(files[1]) ){
        error("Flowcell ID does not match for paired reads of sample ${meta.id} - ${files}")
    }

    // If we cannot read the flowcell ID from the fastq file, then we don't use it
    def sample_lane_id = flowcell ? "${meta.flowcell}.${meta.sample}.${meta.lane}" : "${meta.sample}.${meta.lane}"

    // Don't use a random element for ID, it breaks resuming
    def read_group = "\"@RG\\tID:${sample_lane_id}\\t${CN}PU:${meta.lane}\\tSM:${meta.patient}_${meta.sample}\\tLB:${meta.sample}\\tDS:${params.fasta}\\tPL:${params.seq_platform}\""
    meta  = meta - meta.subMap('lane') + [read_group: read_group.toString()]
    return [ meta, files ]
}

// Parse first line of a FASTQ file, return the flowcell id and lane number.
def flowcellLaneFromFastq(path) {
    // First line of FASTQ file contains sequence identifier plus optional description
    def firstLine = readFirstLineOfFastq(path)
    def flowcell_id = null

    // Expected format from ILLUMINA
    // cf https://en.wikipedia.org/wiki/FASTQ_format#Illumina_sequence_identifiers
    // Five fields:
    // @<instrument>:<lane>:<tile>:<x-pos>:<y-pos>...
    // Seven fields or more (from CASAVA 1.8+):
    // "@<instrument>:<run number>:<flowcell ID>:<lane>:<tile>:<x-pos>:<y-pos>..."

    fields = firstLine ? firstLine.split(':') : []
    if (fields.size() == 5) {
        // Get the instrument name as flowcell ID
        flowcell_id = fields[0].substring(1)
    } else if (fields.size() >= 7) {
        // Get the actual flowcell ID
        flowcell_id = fields[2]
    } else if (fields.size() != 0) {
        log.warn "FASTQ file(${path}): Cannot extract flowcell ID from ${firstLine}"
    }
    return flowcell_id
}

// Get first line of a FASTQ file
def readFirstLineOfFastq(path) {
    def line = null
    try {
        path.withInputStream {
            InputStream gzipStream = new java.util.zip.GZIPInputStream(it)
            Reader decoder = new InputStreamReader(gzipStream, 'ASCII')
            BufferedReader buffered = new BufferedReader(decoder)
            line = buffered.readLine()
            assert line.startsWith('@')
        }
    } catch (Exception e) {
        log.warn "FASTQ file(${path}): Error streaming"
        log.warn "${e.message}"
    }
    return line
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
