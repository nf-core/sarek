//
// This file holds several functions specific to the main.nf workflow in the nf-core/sarek pipeline
//

class WorkflowSarek {

    public static void initialise(params, log, valid_params) {
        genomeExistsError(params, log)

        if (!valid_params['step'].contains(params.step)) {
            log.error "Unknown option: '${params.step}'. Valid options for '--step': ${valid_params['step'].join(', ')}."
            System.exit(1)
        }

        if (params.tools && !valid_params['tools'].contains(params.tools)) {
            log.error "Unknown option: '${params.tools}'. Valid options for '--tools': ${valid_params['tools'].join(', ')}."
            System.exit(1)
        }

        if (params.tools && !valid_params['toolsQcSkip'].contains(params.tools)) {
            log.error "Unknown option: '${params.skip_qc}'. Valid options for '--skip_qc': ${valid_params['toolsQcSkip'].join(', ')}."
            System.exit(1)
        }

        // toolsAnnotate = params.annotate_tools ? params.annotate_tools.split(',').collect{it.trim().toLowerCase().replaceAll('-', '')} : []
        // if (!valid_params['toolsAnnotate'].contains(toolsAnnotate)) {
        //     log.error "Unknown option: '${toolsAnnotate}'. Valid options for '--annotate_tools': ${valid_params['toolsAnnotate'].join(', ')}."
        //     System.exit(1)
        // }

        // if (!valid_params['aligner'].contains(params.aligner)) {
        //     log.error "Unknown option: '${params.aligner}'. Valid options for '--aligner': ${valid_params['aligner'].join(', ')}."
        //     System.exit(1)
        // }
        // if (aligner.contains(',')) exit 1, 'You can choose only one aligner, see --help for more information'

        // // // Check parameters
        // if ((params.ascat_ploidy && !params.ascat_purity) || (!params.ascat_ploidy && params.ascat_purity)) exit 1, 'Please specify both --ascat_purity and --ascat_ploidy, or none of them'
        // if (params.cf_window && params.cf_coeff) exit 1, 'Please specify either --cf_window OR --cf_coeff, but not both of them'
        // if (params.umi && !(params.read_structure1 && params.read_structure2)) exit 1, 'Please specify both --read_structure1 and --read_structure2, when using --umi'
        // if ('mutect2' in tools && !(params.pon)) log.warn "[nf-core/sarek] Mutect2 was requested, but as no panel of normals were given, results will not be optimal"
        // if (params.sentieon) log.warn "[nf-core/sarek] Sentieon will be used, only works if Sentieon is available where nf-core/sarek is run"

        // // Handle input
        // tsv_path = null
        // if (params.input && (hasExtension(params.input, "tsv") || hasExtension(params.input, "vcf") || hasExtension(params.input, "vcf.gz"))) tsv_path = params.input
        // if (params.input && (hasExtension(params.input, "vcf") || hasExtension(params.input, "vcf.gz"))) step = "annotate"

        // if (!params.fasta) {
        //     log.error "Genome fasta file not specified with e.g. '--fasta genome.fa' or via a detectable config file."
        //     System.exit(1)
        // }

        // // If no input file specified, trying to get TSV files corresponding to step in the TSV directory
        // // only for steps preparerecalibration, recalibrate, variantcalling and controlfreec
        // if (!params.input && params.sentieon) {
        //     switch (step) {
        //         case 'mapping': break
        //         case 'recalibrate': tsv_path = "${params.outdir}/preprocessing/tsv/sentieon_deduped.tsv"; break
        //         case 'variantcalling': tsv_path = "${params.outdir}/preprocessing/tsv/sentieon_recalibrated.tsv"; break
        //         case 'annotate': break
        //         default: exit 1, "Unknown step ${step}"
        //     }
        // } else if (!params.input && !params.sentieon && !params.skip_markduplicates) {
        //     switch (step) {
        //         case 'mapping': break
        //         case 'preparerecalibration': tsv_path = "${params.outdir}/preprocessing/tsv/markduplicates_no_table.tsv"; break
        //         case 'recalibrate': tsv_path = "${params.outdir}/preprocessing/tsv/markduplicates.tsv"; break
        //         case 'variantcalling': tsv_path = "${params.outdir}/preprocessing/tsv/recalibrated.tsv"; break
        //         case 'controlfreec': tsv_path = "${params.outdir}/variant_calling/tsv/control-freec_mpileup.tsv"; break
        //         case 'annotate': break
        //         default: exit 1, "Unknown step ${step}"
        //     }
        // } else if (!params.input && !params.sentieon && params.skip_markduplicates) {
        //     switch (step) {
        //         case 'mapping': break
        //         case 'preparerecalibration': tsv_path = "${params.outdir}/preprocessing/tsv/mapped.tsv"; break
        //         case 'recalibrate': tsv_path = "${params.outdir}/preprocessing/tsv/mapped_no_markduplicates.tsv"; break
        //         case 'variantcalling': tsv_path = "${params.outdir}/preprocessing/tsv/recalibrated.tsv"; break
        //         case 'controlfreec': tsv_path = "${params.outdir}/variant_calling/tsv/control-freec_mpileup.tsv"; break
        //         case 'annotate': break
        //         default: exit 1, "Unknown step ${step}"
        //     }
        // }

        // input_sample = Channel.empty()
        // if (tsv_path) {
        //     tsv_file = file(tsv_path)
        //     switch (step) {
        //         case 'mapping': input_sample = extract_fastq(tsv_file); break
        //         case 'preparerecalibration': input_sample = extract_bam(tsv_file); break
        //         case 'recalibrate': input_sample = extract_recal(tsv_file); break
        //         case 'variantcalling': input_sample = extract_bam(tsv_file); break
        //         case 'controlfreec': input_sample = extract_pileup(tsv_file); break
        //         case 'annotate': break
        //         default: exit 1, "Unknown step ${step}"
        //     }
        // } else if (params.input && !hasExtension(params.input, "tsv")) {
        //     log.info "No TSV file"
        //     if (step != 'mapping') exit 1, 'No step other than "mapping" supports a directory as an input'
        //     log.info "Reading ${params.input} directory"
        //     log.warn "[nf-core/sarek] in ${params.input} directory, all fastqs are assuming to be from the same sample, which is assumed to be a germline one"
        //     input_sample = extract_fastq_from_dir(params.input)
        //     tsv_file = params.input  // used in the reports
        // } else if (tsv_path && step == 'annotate') {
        //     log.info "Annotating ${tsv_path}"
        // } else if (step == 'annotate') {
        //     log.info "Trying automatic annotation on files in the VariantCalling/ directory"
        // } else exit 1, 'No sample were defined, see --help'
    }


    // // Check parameter existence
    // public static checkParameterExistence(it, list) {
    //     if (!list.contains(it)) {
    //         log.warn "Unknown parameter: ${it}"
    //         return false
    //     }
    //     return true
    // }

    // // Compare each parameter with a list of parameters
    // public static checkParameterList(list, realList) {
    //     return list.every{ checkParameterExistence(it, realList) }
    // }

    // // Channeling the TSV file containing BAM.
    // // Format is: "patient gender status sample bam bai"
    // public static extractBam(tsvFile) {
    //     Channel.from(tsvFile)
    //         .splitCsv(sep: '\t')
    //         .map { row ->
    //             checkNumberOfItem(row, 6)
    //             def meta = [:]

    //             meta.patient = row[0]
    //             meta.gender  = row[1]
    //             meta.status  = returnStatus(row[2].toInteger())
    //             meta.sample  = row[3]
    //             meta.id      = meta.sample
    //             def bam      = returnFile(row[4])
    //             def bai      = returnFile(row[5])

    //             if (!hasExtension(bam, "bam")) exit 1, "File: ${bam} has the wrong extension. See --help for more information"
    //             if (!hasExtension(bai, "bai")) exit 1, "File: ${bai} has the wrong extension. See --help for more information"

    //             return [meta, bam, bai]
    //         }
    // }

    // // Create a channel of germline FASTQs from a directory pattern: "my_samples/*/"
    // // All FASTQ files in subdirectories are collected and emitted;
    // // they must have _R1_ and _R2_ in their names.
    // // All FASTQ files are assumed to be from the same sample.
    // public static extractFastqFromDir(folder) {
    //     sample = file(folder).getFileName().toString()

    //     fastq = Channel.fromFilePairs(folder + '/*{_R1_,_R2_}*.fastq.gz')
    //         .ifEmpty { error "No directories found matching folder '${folder}'" }

    // // TODO check if flowcellLane_from_fastq is useful or not

    //     fastq = fastq.map{ run, pair ->
    //         def meta = [:]
    //         meta.patient    = sample
    //         meta.sample     = meta.patient
    //         meta.gender     = 'ZZ' // unused
    //         meta.status     = 0    // normal (not tumor)
    //         meta.run        = run
    //         meta.id         = "${meta.sample}-${meta.run}"
    //         def read1       = pair[0]
    //         def read2       = pair[1]
    //         def CN          = params.sequencing_center ? "CN:${params.sequencing_center}\\t" : ""
    //         def read_group  = "\"@RG\\tID:${meta.run}\\t${CN}PU:${meta.run}\\tSM:${meta.sample}\\tLB:${meta.sample}\\tPL:ILLUMINA\""
    //         meta.read_group = read_group

    //         return [meta, [read1, read2]]
    //     }
    // }

    // // Channeling the CSV file containing FASTQ
    // // Format is: "patient,gender,status,sample,lane,fastq1,fastq2"
    // public static extractFastq(csv_file) {
    //     def meta = [:]
    //     csv_file.eachLine{ line ->
    //         def row = line.split(",")
    //         meta.patient    = row[0]
    //         meta.gender     = row[1]
    //         meta.status     = row[2]
    //         meta.sample     = row[3]
    //         meta.run        = row[4]
    //         meta.id         = "${meta.sample}-${meta.run}"
    //         def read1       = file(row[5])
    //         def read2       = file(row[6])
    //         def CN          = params.sequencing_center ? "CN:${params.sequencing_center}\\t" : ''
    //         def read_group  = "\"@RG\\tID:${meta.run}\\t${CN}PU:${meta.run}\\tSM:${meta.sample}\\tLB:${meta.sample}\\tPL:ILLUMINA\""
    //         meta.read_group = read_group
    //         return [meta, [read1, read2]]
    //     }
    // }

    // // // Channeling the TSV file containing mpileup
    // // // Format is: "patient gender status sample pileup"
    // // def extract_pileup(tsvFile) {
    // //     Channel.from(tsvFile)
    // //         .splitCsv(sep: '\t')
    // //         .map { row ->
    // //             checkNumberOfItem(row, 5)
    // //             def idPatient = row[0]
    // //             def gender    = row[1]
    // //             def status    = returnStatus(row[2].toInteger())
    // //             def idSample  = row[3]
    // //             def mpileup   = returnFile(row[4])

    // //             if (!hasExtension(mpileup, "pileup")) exit 1, "File: ${mpileup} has the wrong extension. See --help for more information"

    // //             return [idPatient, gender, status, idSample, mpileup]
    // //         }
    // // }

    // // Channeling the TSV file containing Recalibration Tables.
    // // Format is: "patient gender status sample bam bai recalTable"
    // public static extractRecal(tsvFile) {
    //     Channel.from(tsvFile)
    //         .splitCsv(sep: '\t')
    //         .map { row ->
    //             checkNumberOfItem(row, 7)
    //             def meta = [:]

    //             meta.patient = row[0]
    //             meta.gender  = row[1]
    //             meta.status  = returnStatus(row[2].toInteger())
    //             meta.sample  = row[3]
    //             meta.id      = meta.sample
    //             def bam      = returnFile(row[4])
    //             def bai      = returnFile(row[5])
    //             def table    = returnFile(row[6])

    //             if (!hasExtension(bam, "bam")) exit 1, "File: ${bam} has the wrong extension. See --help for more information"
    //             if (!hasExtension(bai, "bai")) exit 1, "File: ${bai} has the wrong extension. See --help for more information"
    //             if (!hasExtension(table, "recal.table")) exit 1, "File: ${table} has the wrong extension. See --help for more information"

    //             return [meta, bam, bai, table]
    //         }
    // }

    // // // Parse first line of a FASTQ file, return the flowcell id and lane number.
    // // def flowcellLane_from_fastq(path) {
    // //     // expected format:
    // //     // xx:yy:FLOWCELLID:LANE:... (seven fields)
    // //     // or
    // //     // FLOWCELLID:LANE:xx:... (five fields)
    // //     InputStream fileStream = new FileInputStream(path.toFile())
    // //     InputStream gzipStream = new java.util.zip.GZIPInputStream(fileStream)
    // //     Reader decoder = new InputStreamReader(gzipStream, 'ASCII')
    // //     BufferedReader buffered = new BufferedReader(decoder)
    // //     def line = buffered.readLine()
    // //     assert line.startsWith('@')
    // //     line = line.substring(1)
    // //     def fields = line.split(' ')[0].split(':')
    // //     String fcid
    // //     int lane
    // //     if (fields.size() == 7) {
    // //         // CASAVA 1.8+ format
    // //         fcid = fields[2]
    // //         lane = fields[3].toInteger()
    // //     } else if (fields.size() == 5) {
    // //         fcid = fields[0]
    // //         lane = fields[1].toInteger()
    // //     }
    // //     [fcid, lane]
    // // }

    // // Remove .ann .gz and .vcf extension from a VCF file
    // public static reduceVcf(file) {
    //     return file.fileName.toString().minus(".ann").minus(".vcf").minus(".gz")
    // }

    // // Return status [0,1]
    // // 0 == Normal, 1 == Tumor
    // public static returnStatus(it) {
    //     if (!(it in [0, 1])) exit 1, "Status is not recognized in TSV file: ${it}, see --help for more information"
    //     return it
    // }


    //
    // Exit pipeline if incorrect --genome key provided
    //
    private static void genomeExistsError(params, log) {
        if (params.genomes && params.genome && !params.genomes.containsKey(params.genome)) {
            log.error "=============================================================================\n" +
                "  Genome '${params.genome}' not found in any config files provided to the pipeline.\n" +
                "  Currently, the available genome keys are:\n" +
                "  ${params.genomes.keySet().join(", ")}\n" +
                "==================================================================================="
            System.exit(1)
        }
    }
}