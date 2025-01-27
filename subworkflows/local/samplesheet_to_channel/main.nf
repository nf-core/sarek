workflow  SAMPLESHEET_TO_CHANNEL{

    take:
    ch_from_samplesheet             //
    aligner                         //
    ascat_alleles                   //
    ascat_loci                      //
    ascat_loci_gc                   //
    ascat_loci_rt                   //
    bcftools_annotations            //
    bcftools_annotations_tbi        //
    bcftools_header_lines           //
    build_only_index                //
    dbsnp                           //
    fasta                           //
    germline_resource               //
    intervals                       //
    joint_germline                  //
    joint_mutect2                   //
    known_indels                    //
    known_snps                      //
    no_intervals                    //
    pon                             //
    sentieon_dnascope_emit_mode     //
    sentieon_haplotyper_emit_mode   //
    seq_center                      //
    seq_platform                    //
    skip_tools                      //
    snpeff_cache                    //
    snpeff_db                       //
    step                            //
    tools                           //
    umi_read_structure              //
    wes                             //

    main:
    ch_from_samplesheet.dump(tag:"ch_from_samplesheet")
    input_sample = ch_from_samplesheet.map{ meta, fastq_1, fastq_2, spring_1, spring_2, table, cram, crai, bam, bai, vcf, variantcaller ->
        // generate patient_sample key to group lanes together
        [ meta.patient + meta.sample, [meta, fastq_1, fastq_2, spring_1, spring_2, table, cram, crai, bam, bai, vcf, variantcaller] ]
    }.tap{ ch_with_patient_sample } // save the channel
    .groupTuple() //group by patient_sample to get all lanes
    .map { patient_sample, ch_items ->
        // get number of lanes per sample
        [ patient_sample, ch_items.size() ]
    }.combine(ch_with_patient_sample, by: 0) // for each entry add numLanes
    .map { patient_sample, num_lanes, ch_items ->
        (meta, fastq_1, fastq_2, spring_1, spring_2, table, cram, crai, bam, bai, vcf, variantcaller) = ch_items
        if (meta.lane && fastq_2) {
            meta = meta + [id: "${meta.sample}-${meta.lane}".toString(), data_type: "fastq_gz", num_lanes: num_lanes.toInteger(), size: 1]

            if (step == 'mapping') return [ meta, [ fastq_1, fastq_2 ] ]
            else {
                error("Samplesheet contains fastq files but step is `$step`. Please check your samplesheet or adjust the step parameter.\nhttps://nf-co.re/sarek/usage#input-samplesheet-configurations")
            }

        // start from TWO spring-files - one with R1 and one with R2
        } else if (meta.lane && spring_1 && spring_2) {
            meta = meta + [id: "${meta.sample}-${meta.lane}".toString(), data_type: "two_fastq_gz_spring", num_lanes: num_lanes.toInteger(), size: 1]

            if (step == 'mapping') return [ meta, [ spring_1, spring_2 ] ]
            else {
                error("Samplesheet contains spring files (in columns `spring_1` and `spring_2`) but step is `$step`. Please check your samplesheet or adjust the step parameter.\nhttps://nf-co.re/sarek/usage#input-samplesheet-configurations")
            }

        // start from ONE spring-file containing both R1 and R2
        } else if (meta.lane && spring_1 && !spring_2) {
            meta = meta + [id: "${meta.sample}-${meta.lane}".toString(), data_type: "one_fastq_gz_spring", num_lanes: num_lanes.toInteger(), size: 1]

            if (step == 'mapping') return [ meta, [ spring_1 ] ]
            else {
                error("Samplesheet contains a spring file (in columns `spring_1`) but step is `$step`. Please check your samplesheet or adjust the step parameter.\nhttps://nf-co.re/sarek/usage#input-samplesheet-configurations")
            }

        // start from BAM
        } else if (meta.lane && bam) {
            if (step != 'mapping' && !bai) {
                error("BAM index (bai) should be provided.")
            }
            meta            = meta + [id: "${meta.sample}-${meta.lane}".toString()]
            def CN          = seq_center ? "CN:${seq_center}\\t" : ''
            def read_group  = "\"@RG\\tID:${meta.sample}_${meta.lane}\\t${CN}PU:${meta.lane}\\tSM:${meta.patient}_${meta.sample}\\tLB:${meta.sample}\\tDS:${fasta}\\tPL:${seq_platform}\""

            meta            = meta - meta.subMap('lane') + [num_lanes: num_lanes.toInteger(), read_group: read_group.toString(), data_type: 'bam', size: 1]

            if (step != 'annotate') return [ meta - meta.subMap('lane'), bam, bai ]
            else {
                error("Samplesheet contains bam files but step is `annotate`. The pipeline is expecting vcf files for the annotation. Please check your samplesheet or adjust the step parameter.\nhttps://nf-co.re/sarek/usage#input-samplesheet-configurations")
            }

        // recalibration
        } else if (table && cram) {
            meta = meta + [id: meta.sample, data_type: 'cram']

            if (!(step == 'mapping' || step == 'annotate')) return [ meta - meta.subMap('lane'), cram, crai, table ]
            else {
                error("Samplesheet contains cram files but step is `$step`. Please check your samplesheet or adjust the step parameter.\nhttps://nf-co.re/sarek/usage#input-samplesheet-configurations")
            }

        // recalibration when skipping MarkDuplicates
        } else if (table && bam) {
            meta = meta + [id: meta.sample, data_type: 'bam']

            if (!(step == 'mapping' || step == 'annotate')) return [ meta - meta.subMap('lane'), bam, bai, table ]
            else {
                error("Samplesheet contains bam files but step is `$step`. Please check your samplesheet or adjust the step parameter.\nhttps://nf-co.re/sarek/usage#input-samplesheet-configurations")
            }

        // prepare_recalibration or variant_calling
        } else if (cram) {
            meta = meta + [id: meta.sample, data_type: 'cram']

            if (!(step == 'mapping' || step == 'annotate')) return [ meta - meta.subMap('lane'), cram, crai ]
            else {
                error("Samplesheet contains cram files but step is `$step`. Please check your samplesheet or adjust the step parameter.\nhttps://nf-co.re/sarek/usage#input-samplesheet-configurations")
            }

        // prepare_recalibration when skipping MarkDuplicates or `--step markduplicates`
        } else if (bam) {
            meta = meta + [id: meta.sample, data_type: 'bam']

            if (!(step == 'mapping' || step == 'annotate')) return [ meta - meta.subMap('lane'), bam, bai ]
            else {
                error("Samplesheet contains bam files but step is `$step`. Please check your samplesheet or adjust the step parameter.\nhttps://nf-co.re/sarek/usage#input-samplesheet-configurations")
            }

        // annotation
        } else if (vcf) {
            meta = meta + [id: meta.sample, data_type: 'vcf', variantcaller: variantcaller ?: '']

            if (step == 'annotate') return [ meta - meta.subMap('lane'), vcf ]
            else {
                error("Samplesheet contains vcf files but step is `$step`. Please check your samplesheet or adjust the step parameter.\nhttps://nf-co.re/sarek/usage#input-samplesheet-configurations")
            }
        } else {
            error("Missing or unknown field in csv file header. Please check your samplesheet")
        }
    }

    if (step != 'annotate' && tools && !build_only_index) {
    // Two checks for ensuring that the pipeline stops with a meaningful error message if
    // 1. the sample-sheet only contains normal-samples, but some of the requested tools require tumor-samples, and
    // 2. the sample-sheet only contains tumor-samples, but some of the requested tools require normal-samples.
    input_sample.filter{ it[0].status == 1 }.ifEmpty{ // In this case, the sample-sheet contains no tumor-samples
        if (!build_only_index) {
            def tools_tumor = ['ascat', 'controlfreec', 'mutect2', 'msisensorpro']
            def tools_tumor_asked = []
            tools_tumor.each{ tool ->
                if (tools.split(',').contains(tool)) tools_tumor_asked.add(tool)
            }
            if (!tools_tumor_asked.isEmpty()) {
                error('The sample-sheet only contains normal-samples, but the following tools, which were requested with "--tools", expect at least one tumor-sample : ' + tools_tumor_asked.join(", "))
            }
        }
    }

    input_sample.filter{ it[0].status == 0 }.ifEmpty{ // In this case, the sample-sheet contains no normal/germline-samples
        def tools_requiring_normal_samples = ['ascat', 'deepvariant', 'haplotypecaller', 'msisensorpro']
        def requested_tools_requiring_normal_samples = []
        tools_requiring_normal_samples.each{ tool_requiring_normal_samples ->
            if (tools.split(',').contains(tool_requiring_normal_samples)) requested_tools_requiring_normal_samples.add(tool_requiring_normal_samples)
        }
        if (!requested_tools_requiring_normal_samples.isEmpty()) {
            error('The sample-sheet only contains tumor-samples, but the following tools, which were requested by the option "tools", expect at least one normal-sample : ' + requested_tools_requiring_normal_samples.join(", "))
        }
    }
    }

    // Fails when wrongfull extension for intervals file
    if (wes && !step == 'annotate') {
        if (intervals && !intervals.endsWith("bed"))  error("Target file specified with `--intervals` must be in BED format for targeted data")
        else log.warn("Intervals file was provided without parameter `--wes`: Pipeline will assume this is Whole-Genome-Sequencing data.")
    } else if (intervals && !intervals.endsWith("bed") && !intervals.endsWith("list")) error("Intervals file must end with .bed, .list, or .interval_list")

    if (step == 'mapping' && aligner.contains("dragmap") && !(skip_tools && skip_tools.split(',').contains("baserecalibrator"))) {
        log.warn("DragMap was specified as aligner. Base recalibration is not contained in --skip_tools. It is recommended to skip baserecalibration when using DragMap\nhttps://gatk.broadinstitute.org/hc/en-us/articles/4407897446939--How-to-Run-germline-single-sample-short-variant-discovery-in-DRAGEN-mode")
    }

    if (step == 'mapping' && aligner.contains("sentieon-bwamem") && umi_read_structure) {
        error("Sentieon BWA is currently not compatible with FGBio UMI handeling. Please choose a different aligner.")
    }

    if (tools && tools.split(',').contains("sentieon_haplotyper") && joint_germline && (!sentieon_haplotyper_emit_mode || !(sentieon_haplotyper_emit_mode.contains('gvcf')))) {
        error("When setting the option `--joint_germline` and including `sentieon_haplotyper` among the requested tools, please set `--sentieon_haplotyper_emit_mode` to include `gvcf`.")
    }

    // Fails or warns when missing files or params for ascat
    if (tools && tools.split(',').contains('ascat')) {
        if (!ascat_alleles) {
            error("No allele files were provided for running ASCAT. Please provide a zip folder with allele files.")
        }
        if (!ascat_loci) {
            error("No loci files were provided for running ASCAT. Please provide a zip folder with loci files.")
        }
        if (!ascat_loci_gc && !ascat_loci_rt) {
            log.warn("No LogRCorrection performed in ASCAT. For LogRCorrection to run, please provide either loci gc files or both loci gc files and loci rt files.")
        }
        if (wes) {
            log.warn("Default reference files not suited for running ASCAT on WES data. It's recommended to use the reference files provided here: https://github.com/Wedge-lab/battenberg#required-reference-files")
        }
    }

    // Warns when missing files or params for mutect2
    if (tools && tools.split(',').contains('mutect2')) {
        if (!pon) {
            log.warn("No Panel-of-normal was specified for Mutect2.\nIt is highly recommended to use one: https://gatk.broadinstitute.org/hc/en-us/articles/5358911630107-Mutect2\nFor more information on how to create one: https://gatk.broadinstitute.org/hc/en-us/articles/5358921041947-CreateSomaticPanelOfNormals-BETA-")
        }
        if (!germline_resource) {
            log.warn("If Mutect2 is specified without a germline resource, no filtering will be done.\nIt is recommended to use one: https://gatk.broadinstitute.org/hc/en-us/articles/5358911630107-Mutect2")
        }
        if (pon && pon.contains("/Homo_sapiens/GATK/GRCh38/Annotation/GATKBundle/1000g_pon.hg38.vcf.gz")) {
            log.warn("The default Panel-of-Normals provided by GATK is used for Mutect2.\nIt is highly recommended to generate one from normal samples that are technical similar to the tumor ones.\nFor more information: https://gatk.broadinstitute.org/hc/en-us/articles/360035890631-Panel-of-Normals-PON-")
        }
    }

    // Fails when missing resources for baserecalibrator
    // Warns when missing resources for haplotypecaller
    if (!dbsnp && !known_indels) {
        if (step in ['mapping', 'markduplicates', 'prepare_recalibration', 'recalibrate'] && (!skip_tools || (skip_tools && !skip_tools.split(',').contains('baserecalibrator')))) {
            error("Base quality score recalibration requires at least one resource file. Please provide at least one of `--dbsnp` or `--known_indels`\nYou can skip this step in the workflow by adding `--skip_tools baserecalibrator` to the command.")
        }
        if (tools && (tools.split(',').contains('haplotypecaller') || tools.split(',').contains('sentieon_haplotyper') || tools.split(',').contains('sentieon_dnascope'))) {
            log.warn "If GATK's Haplotypecaller, Sentieon's Dnascope or Sentieon's Haplotyper is specified, without `--dbsnp` or `--known_indels no filtering will be done. For filtering, please provide at least one of `--dbsnp` or `--known_indels`.\nFor more information see FilterVariantTranches (single-sample, default): https://gatk.broadinstitute.org/hc/en-us/articles/5358928898971-FilterVariantTranches\nFor more information see VariantRecalibration (--joint_germline): https://gatk.broadinstitute.org/hc/en-us/articles/5358906115227-VariantRecalibrator\nFor more information on GATK Best practice germline variant calling: https://gatk.broadinstitute.org/hc/en-us/articles/360035535932-Germline-short-variant-discovery-SNPs-Indels-"
        }
    }
    if (joint_germline && (!tools || !(tools.split(',').contains('haplotypecaller') || tools.split(',').contains('sentieon_haplotyper') || tools.split(',').contains('sentieon_dnascope')))) {
        error("The GATK's Haplotypecaller, Sentieon's Dnascope or Sentieon's Haplotyper should be specified as one of the tools when doing joint germline variant calling.) ")
    }

    if (
        tools &&
        (
            tools.split(',').contains('haplotypecaller') ||
            tools.split(',').contains('sentieon_haplotyper') ||
            tools.split(',').contains('sentieon_dnascope')
        ) &&
        joint_germline &&
            ( !dbsnp || !known_indels || !known_snps || no_intervals )
    ) {
        log.warn("""If GATK's Haplotypecaller, Sentieon's Dnascope and/or Sentieon's Haplotyper is specified, but without `--dbsnp`, `--known_snps`, `--known_indels` or the associated resource labels (ie `known_snps_vqsr`), no variant recalibration will be done. For recalibration you must provide all of these resources.\nFor more information see VariantRecalibration: https://gatk.broadinstitute.org/hc/en-us/articles/5358906115227-VariantRecalibrator \n\
Joint germline variant calling also requires intervals in order to genotype the samples. As a result, if `--no_intervals` is set to `true` the joint germline variant calling will not be performed.""")
    }

    if (tools &&
        tools.split(',').contains('sentieon_dnascope') && joint_germline &&
        ( !sentieon_dnascope_emit_mode || !sentieon_dnascope_emit_mode.split(',').contains('gvcf') )
    ) {
        error("When using Sentieon Dnascope for joint-germline variant-calling the option `--sentieon_dnascope_emit_mode` has to include `gvcf`.")
    }

    if (tools &&
        tools.split(',').contains('sentieon_haplotyper') && joint_germline &&
        ( !sentieon_haplotyper_emit_mode || !sentieon_haplotyper_emit_mode.split(',').contains('gvcf') )
    ) {
        error("When using Sentieon Haplotyper for joint-germline variant-calling the option `--sentieon_haplotyper_emit_mode` has to include `gvcf`.")
    }


    // Fails when --joint_mutect2 is used without enabling mutect2
    if (joint_mutect2 && (!tools || !tools.split(',').contains('mutect2'))) {
        error("The mutect2 should be specified as one of the tools when doing joint somatic variant calling with Mutect2. (The mutect2 could be specified by adding `--tools mutect2` to the nextflow command.)")
    }

    // Fails when missing tools for variant_calling or annotate
    if ((step == 'variant_calling' || step == 'annotate') && !tools) {
        error("Please specify at least one tool when using `--step ${step}`.\nhttps://nf-co.re/sarek/parameters#tools")
    }

    // Fails when missing sex information for CNV tools
    if (tools && (tools.split(',').contains('ascat') || tools.split(',').contains('controlfreec'))) {
        input_sample.map{
            if (it[0].sex == 'NA' ) {
                error("Please specify sex information for each sample in your samplesheet when using '--tools' with 'ascat' or 'controlfreec'.\nhttps://nf-co.re/sarek/usage#input-samplesheet-configurations")
            }
        }
    }

    // Fails when bcftools annotate is used but no files are supplied
    if (tools && tools.split(',').contains('bcfann') && !(bcftools_annotations && bcftools_annotations_tbi && bcftools_header_lines)) {
        error("Please specify --bcftools_annotations, --bcftools_annotations_tbi, and --bcftools_header_lines, when using BCFTools annotations")
    }

    // Fails when snpeff annotation is enabled but snpeff_db is not specified
    if ((snpeff_cache && tools && (tools.split(',').contains("snpeff") || tools.split(',').contains('merge'))) &&
        !snpeff_db) {
        error("Please specify --snpeff_db")
    }

    emit:
    input_sample
    }
