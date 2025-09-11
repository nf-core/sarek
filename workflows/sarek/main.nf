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
include { CHANNEL_VARIANT_CALLING_CREATE_CSV                } from '../../subworkflows/local/channel_variant_calling_create_csv'

// Convert BAM files to FASTQ files
include { BAM_CONVERT_SAMTOOLS as CONVERT_FASTQ_INPUT       } from '../../subworkflows/local/bam_convert_samtools'

// Convert fastq.gz.spring files to fastq.gz files
include { SPRING_DECOMPRESS as SPRING_DECOMPRESS_TO_R1_FQ   } from '../../modules/nf-core/spring/decompress'
include { SPRING_DECOMPRESS as SPRING_DECOMPRESS_TO_R2_FQ   } from '../../modules/nf-core/spring/decompress'
include { SPRING_DECOMPRESS as SPRING_DECOMPRESS_TO_FQ_PAIR } from '../../modules/nf-core/spring/decompress'

// Run FASTQC
include { FASTQC                                            } from '../../modules/nf-core/fastqc'

// QC on CRAM
include { CRAM_SAMPLEQC                                     } from '../../subworkflows/local/cram_sampleqc'

// Preprocessing
include { FASTQ_PREPROCESS_GATK                             } from '../../subworkflows/local/fastq_preprocess_gatk'
include { FASTQ_PREPROCESS_PARABRICKS                       } from '../../subworkflows/local/fastq_preprocess_parabricks'

// Variant calling on a single normal sample
include { BAM_VARIANT_CALLING_GERMLINE_ALL                  } from '../../subworkflows/local/bam_variant_calling_germline_all'

// Variant calling on a single tumor sample
include { BAM_VARIANT_CALLING_TUMOR_ONLY_ALL                } from '../../subworkflows/local/bam_variant_calling_tumor_only_all'

// Variant calling on tumor/normal pair
include { BAM_VARIANT_CALLING_SOMATIC_ALL                   } from '../../subworkflows/local/bam_variant_calling_somatic_all'

// POST VARIANTCALLING: e.g. merging
include { POST_VARIANTCALLING                               } from '../../subworkflows/local/post_variantcalling'

// QC on VCF files
include { VCF_QC_BCFTOOLS_VCFTOOLS                          } from '../../subworkflows/local/vcf_qc_bcftools_vcftools'

// Annotation
include { VCF_ANNOTATE_ALL                                  } from '../../subworkflows/local/vcf_annotate_all'

// MULTIQC
include { MULTIQC                                           } from '../../modules/nf-core/multiqc'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow SAREK {
    take:
        input_sample
        allele_files
        aligner
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
        versions

    main:

    // To gather all QC reports for MultiQC
    ch_multiqc_files = Channel.empty()
    multiqc_report   = Channel.empty()
    reports          = Channel.empty()

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

        versions = versions.mix(SPRING_DECOMPRESS_TO_R1_FQ.out.versions)
        versions = versions.mix(SPRING_DECOMPRESS_TO_R2_FQ.out.versions)
        versions = versions.mix(SPRING_DECOMPRESS_TO_FQ_PAIR.out.versions)

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

        versions = versions.mix(CONVERT_FASTQ_INPUT.out.versions)

        // Gather fastq (inputed or converted)
        // Theorically this could work on mixed input (fastq for one sample and bam for another)
        // But not sure how to handle that with the samplesheet
        // Or if we really want users to be able to do that
        input_fastq = fastq_gz.mix(CONVERT_FASTQ_INPUT.out.reads).mix(one_fastq_gz_from_spring).mix(two_fastq_gz_from_spring)

        // QC
        // `--skip_tools fastqc` to skip fastqc
        if (!(params.skip_tools && params.skip_tools.split(',').contains('fastqc'))) {
            FASTQC(input_fastq)

            reports = reports.mix(FASTQC.out.zip.collect{ _meta, logs -> logs })
            versions = versions.mix(FASTQC.out.versions)
        }
    }
    else {
        input_fastq = Channel.empty().mix( input_sample )
    }

    if (params.step in ['mapping', 'markduplicates', 'prepare_recalibration', 'recalibrate']) {

        if (aligner == 'parabricks') {
            // PREPROCESSING WITH PARABRICKS
            FASTQ_PREPROCESS_PARABRICKS(
                input_fastq,
                fasta,
                index_alignment,
                intervals_and_num_intervals,
                known_sites_indels,
                Channel.value("cram")
            )

            // Gather preprocessing output
            cram_variant_calling = Channel.empty()
            cram_variant_calling = cram_variant_calling.mix(FASTQ_PREPROCESS_PARABRICKS.out.cram)

            // Gather used softwares versions
            reports = reports.mix(FASTQ_PREPROCESS_PARABRICKS.out.reports)
            versions = versions.mix(FASTQ_PREPROCESS_PARABRICKS.out.versions)
        } else {
            // PREPROCESSING
            FASTQ_PREPROCESS_GATK(
                input_fastq,
                input_sample,
                dict,
                fasta,
                fasta_fai,
                index_alignment,
                intervals_and_num_intervals,
                intervals_for_preprocessing,
                known_sites_indels,
                known_sites_indels_tbi)

            // Gather preprocessing output
            cram_variant_calling = Channel.empty()
            cram_variant_calling = cram_variant_calling.mix(FASTQ_PREPROCESS_GATK.out.cram_variant_calling)

            // Gather used softwares versions
            reports = reports.mix(FASTQ_PREPROCESS_GATK.out.reports)
            versions = versions.mix(FASTQ_PREPROCESS_GATK.out.versions)
        }

    }

    if (params.step == 'variant_calling') {

        cram_variant_calling = Channel.empty().mix( input_sample )

    }

    if (params.step == 'annotate') {

        cram_variant_calling = Channel.empty()

    }

    // RUN CRAM QC on the recalibrated CRAM files or when starting from step variant calling. NGSCheckmate should be run also on non-recalibrated CRAM files
    CRAM_SAMPLEQC(cram_variant_calling,
        ngscheckmate_bed,
        fasta,
        params.skip_tools && params.skip_tools.split(',').contains('baserecalibrator'),
        intervals_for_preprocessing)

    reports  = reports.mix(CRAM_SAMPLEQC.out.reports)
    versions = versions.mix(CRAM_SAMPLEQC.out.versions)

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

                meta.id            = "${tumor[1].sample}_vs_${normal[1].sample}".toString()
                meta.normal_id     = normal[1].sample
                meta.patient       = normal[0]
                meta.sex           = normal[1].sex
                meta.tumor_id      = tumor[1].sample
                meta.contamination = tumor[1].contamination ?: 0.5

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
        POST_VARIANTCALLING(params.tools,
                cram_variant_calling_status_normal,
                BAM_VARIANT_CALLING_GERMLINE_ALL.out.vcf_all,
                cram_variant_calling_tumor_only,
                BAM_VARIANT_CALLING_TUMOR_ONLY_ALL.out.vcf_all,
                cram_variant_calling_pair,
                BAM_VARIANT_CALLING_SOMATIC_ALL.out.vcf_all,
                fasta,
                fasta_fai,
                params.concatenate_vcfs,
                params.normalize_vcfs,
                params.varlociraptor_chunk_size,
            )

        // Gather vcf files for annotation and QC
        vcf_to_annotate = Channel.empty()

        // Check if normalization is requested
        if (params.normalize_vcfs) {
            vcf_to_annotate = vcf_to_annotate.mix(POST_VARIANTCALLING.out.vcfs)
        } else {
            // If not normalized, gather existing VCFs
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
        }

        // QC
        VCF_QC_BCFTOOLS_VCFTOOLS(vcf_to_annotate, intervals_bed_combined)

        reports = reports.mix(VCF_QC_BCFTOOLS_VCFTOOLS.out.bcftools_stats.collect{ _meta, stats -> [ stats ] })
        reports = reports.mix(VCF_QC_BCFTOOLS_VCFTOOLS.out.vcftools_tstv_counts.collect{ _meta, counts -> [ counts ] })
        reports = reports.mix(VCF_QC_BCFTOOLS_VCFTOOLS.out.vcftools_tstv_qual.collect{ _meta, qual -> [ qual ] })
        reports = reports.mix(VCF_QC_BCFTOOLS_VCFTOOLS.out.vcftools_filter_summary.collect{ _meta, summary -> [ summary ] })
        reports = reports.mix(BAM_VARIANT_CALLING_GERMLINE_ALL.out.out_indexcov.collect{ _meta, indexcov -> indexcov.flatten() })
        reports = reports.mix(BAM_VARIANT_CALLING_SOMATIC_ALL.out.out_indexcov.collect{ _meta, indexcov -> indexcov.flatten() })

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
        .collectFile(storeDir: "${params.outdir}/pipeline_info", name: 'nf_core_'  +  'sarek_software_'  + 'mqc_'  + 'versions.yml', sort: true, newLine: true)
    }

    //
    // MODULE: MultiQC
    //
    if (!(params.skip_tools && params.skip_tools.split(',').contains('multiqc'))) {
        ch_multiqc_config = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
        ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config, checkIfExists: true) : Channel.empty()
        ch_multiqc_logo = params.multiqc_logo ? Channel.fromPath(params.multiqc_logo, checkIfExists: true) : Channel.empty()

        summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
        ch_workflow_summary = Channel.value(paramsSummaryMultiqc(summary_params))
        ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
        ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
        ch_methods_description = Channel.value(methodsDescriptionText(ch_multiqc_custom_methods_description))

        ch_multiqc_files = ch_multiqc_files.mix(version_yaml)
        ch_multiqc_files = ch_multiqc_files.mix(reports)
        ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml', sort: true))

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
    def sample_lane_id = flowcell ? "${flowcell}.${meta.sample}.${meta.lane}" : "${meta.sample}.${meta.lane}"

    // Don't use a random element for ID, it breaks resuming
    def read_group = params.umi_read_structure ?
        "\"@RG\\tID:${meta.sample}\\t${CN}PU:consensus\\tSM:${meta.patient}_${meta.sample}\\tLB:${meta.sample}\\tDS:${params.fasta}\\tPL:${params.seq_platform}\"" :
        "\"@RG\\tID:${sample_lane_id}\\t${CN}PU:${meta.lane}\\tSM:${meta.patient}_${meta.sample}\\tLB:${meta.sample}\\tDS:${params.fasta}\\tPL:${params.seq_platform}\""

    meta  = meta - meta.subMap('lane') + [read_group: read_group.toString(), sample_lane_id: sample_lane_id.toString()]
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

    def fields = firstLine ? firstLine.split(':') : []
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
