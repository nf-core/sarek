/*
 * This file holds several functions used to perform operation in Sarek
 */
 
// Check if a row has the expected number of item
def check_number_of_item(row, number) {
    if (row.size() != number) exit 1, "Malformed row in TSV file: ${row}, see --help for more information"
    return true
}

// Check parameter existence
def check_parameter_existence(it, list) {
    if (!list.contains(it)) {
        log.warn "Unknown parameter: ${it}"
        return false
    }
    return true
}

// Compare each parameter with a list of parameters
def check_parameter_list(list, realList) {
    return list.every{ check_parameter_existence(it, realList) }
}

// Define list of available tools to annotate
def define_anno_list() {
    return [
        'haplotypecaller',
        'manta',
        'mutect2',
        'strelka',
        'tiddit'
    ]
}

// Define list of skipable QC tools
def define_skip_qc_list() {
    return [
        'bamqc',
        'baserecalibrator',
        'bcftools',
        'documentation',
        'fastqc',
        'markduplicates',
        'multiqc',
        'samtools',
        'sentieon',
        'vcftools',
        'versions'
    ]
}

// Define list of available step
def define_step_list() {
    return [
        'annotate',
        'controlfreec',
        'mapping',
        'preparerecalibration',
        'recalibrate',
        'variantcalling'
    ]
}

// Define list of available tools
def define_tool_list() {
    return [
        'ascat',
        'cnvkit',
        'controlfreec',
        'dnascope',
        'dnaseq',
        'freebayes',
        'haplotypecaller',
        'manta',
        'merge',
        'mpileup',
        'mutect2',
        'snpeff',
        'strelka',
        'tiddit',
        'tnscope',
        'vep',
        'msisensor'
    ]
}

// Channeling the TSV file containing BAM.
// Format is: "patient gender status sample bam bai"
def extract_bam(tsvFile) {
    Channel.from(tsvFile)
        .splitCsv(sep: '\t')
        .map { row ->
            check_number_of_item(row, 6)
            def meta = [:]

            meta.patient = row[0]
            meta.gender  = row[1]
            meta.status  = return_status(row[2].toInteger())
            meta.sample  = row[3]
            meta.id      = meta.sample
            def bam      = return_file(row[4])
            def bai      = return_file(row[5])

            if (!has_extension(bam, "bam")) exit 1, "File: ${bam} has the wrong extension. See --help for more information"
            if (!has_extension(bai, "bai")) exit 1, "File: ${bai} has the wrong extension. See --help for more information"

            return [meta, bam, bai]
        }
}

// Create a channel of germline FASTQs from a directory pattern: "my_samples/*/"
// All FASTQ files in subdirectories are collected and emitted;
// they must have _R1_ and _R2_ in their names.
// All FASTQ files are assumed to be from the same sample.
def extract_fastq_from_dir(folder) {
    sample = file(folder).getFileName().toString()

    fastq = Channel.fromFilePairs(folder + '/*{_R1_,_R2_}*.fastq.gz')
        .ifEmpty { error "No directories found matching folder '${folder}'" }

// TODO check if flowcellLane_from_fastq is useful or not

    fastq = fastq.map{ run, pair ->
        def meta = [:]
        meta.patient    = sample
        meta.sample     = meta.patient
        meta.gender     = 'ZZ' // unused
        meta.status     = 0    // normal (not tumor)
        meta.run        = run
        meta.id         = "${meta.sample}-${meta.run}"
        def read1       = pair[0]
        def read2       = pair[1]
        def CN          = params.sequencing_center ? "CN:${params.sequencing_center}\\t" : ""
        def read_group  = "\"@RG\\tID:${meta.run}\\t${CN}PU:${meta.run}\\tSM:${meta.sample}\\tLB:${meta.sample}\\tPL:ILLUMINA\""
        meta.read_group = read_group

        return [meta, [read1, read2]]
    }
}

// Channeling the TSV file containing FASTQ or BAM
// Format is: "patient gender status sample lane fastq1 fastq2"
// or: "patient gender status sample lane bam"
def extract_fastq(tsvFile) {
    Channel.from(tsvFile)
        .splitCsv(sep: '\t')
        .map { row ->
            def meta = [:]
            meta.patient    = row[0]
            meta.gender     = row[1]
            meta.status     = return_status(row[2].toInteger())
            meta.sample     = row[3]
            meta.run        = row[4]
            meta.id         = "${meta.sample}-${meta.run}"
            def read1       = return_file(row[5])
            def read2       = "null"
            def CN          = params.sequencing_center ? "CN:${params.sequencing_center}\\t" : ""
            def read_group  = "\"@RG\\tID:${meta.run}\\t${CN}PU:${meta.run}\\tSM:${meta.sample}\\tLB:${meta.sample}\\tPL:ILLUMINA\""
            meta.read_group = read_group

            if (has_extension(read1, "fastq.gz") || has_extension(read1, "fq.gz") || has_extension(read1, "fastq") || has_extension(read1, "fq")) {
                check_number_of_item(row, 7)
                read2 = return_file(row[6])
            if (!has_extension(read2, "fastq.gz") && !has_extension(read2, "fq.gz")  && !has_extension(read2, "fastq") && !has_extension(read2, "fq")) exit 1, "File: ${file2} has the wrong extension. See --help for more information"
            if (has_extension(read1, "fastq") || has_extension(read1, "fq") || has_extension(read2, "fastq") || has_extension(read2, "fq")) {
                exit 1, "We do recommend to use gziped fastq file to help you reduce your data footprint."
            }
        }
        else if (has_extension(read1, "bam")) check_number_of_item(row, 6)
        else exit 1, "No recognisable extention for input file: ${read1}"

        return [meta, [read1, read2]]
    }
}

// // Channeling the TSV file containing mpileup
// // Format is: "patient gender status sample pileup"
// def extract_pileup(tsvFile) {
//     Channel.from(tsvFile)
//         .splitCsv(sep: '\t')
//         .map { row ->
//             check_number_of_item(row, 5)
//             def idPatient = row[0]
//             def gender    = row[1]
//             def status    = return_status(row[2].toInteger())
//             def idSample  = row[3]
//             def mpileup   = return_file(row[4])

//             if (!has_extension(mpileup, "pileup")) exit 1, "File: ${mpileup} has the wrong extension. See --help for more information"

//             return [idPatient, gender, status, idSample, mpileup]
//         }
// }

// Channeling the TSV file containing Recalibration Tables.
// Format is: "patient gender status sample bam bai recalTable"
def extract_recal(tsvFile) {
    Channel.from(tsvFile)
        .splitCsv(sep: '\t')
        .map { row ->
            check_number_of_item(row, 7)
            def meta = [:]

            meta.patient = row[0]
            meta.gender  = row[1]
            meta.status  = return_status(row[2].toInteger())
            meta.sample  = row[3]
            meta.id      = meta.sample
            def bam      = return_file(row[4])
            def bai      = return_file(row[5])
            def table    = return_file(row[6])

            if (!has_extension(bam, "bam")) exit 1, "File: ${bam} has the wrong extension. See --help for more information"
            if (!has_extension(bai, "bai")) exit 1, "File: ${bai} has the wrong extension. See --help for more information"
            if (!has_extension(table, "recal.table")) exit 1, "File: ${table} has the wrong extension. See --help for more information"

            return [meta, bam, bai, table]
        }
}

// // Parse first line of a FASTQ file, return the flowcell id and lane number.
// def flowcellLane_from_fastq(path) {
//     // expected format:
//     // xx:yy:FLOWCELLID:LANE:... (seven fields)
//     // or
//     // FLOWCELLID:LANE:xx:... (five fields)
//     InputStream fileStream = new FileInputStream(path.toFile())
//     InputStream gzipStream = new java.util.zip.GZIPInputStream(fileStream)
//     Reader decoder = new InputStreamReader(gzipStream, 'ASCII')
//     BufferedReader buffered = new BufferedReader(decoder)
//     def line = buffered.readLine()
//     assert line.startsWith('@')
//     line = line.substring(1)
//     def fields = line.split(' ')[0].split(':')
//     String fcid
//     int lane
//     if (fields.size() == 7) {
//         // CASAVA 1.8+ format
//         fcid = fields[2]
//         lane = fields[3].toInteger()
//     } else if (fields.size() == 5) {
//         fcid = fields[0]
//         lane = fields[1].toInteger()
//     }
//     [fcid, lane]
// }

// Check file extension
def has_extension(it, extension) {
    it.toString().toLowerCase().endsWith(extension.toLowerCase())
}

// Return file if it exists
def return_file(it) {
    if (!file(it).exists()) exit 1, "Missing file in TSV file: ${it}, see --help for more information"
    return file(it)
}

// Remove .ann .gz and .vcf extension from a VCF file
def reduce_vcf(file) {
    return file.fileName.toString().minus(".ann").minus(".vcf").minus(".gz")
}

// Return status [0,1]
// 0 == Normal, 1 == Tumor
def return_status(it) {
    if (!(it in [0, 1])) exit 1, "Status is not recognized in TSV file: ${it}, see --help for more information"
    return it
}

/*
 * nf-core core functions
 */

/*
 * Extract name of software tool from process name using $task.process
 */
def getSoftwareName(task_process) {
    return task_process.tokenize(':')[-1].tokenize('_')[0].toLowerCase()
}

/*
 * Function to initialise default values and to generate a Groovy Map of available options for nf-core modules
 */
def initOptions(Map args) {
    def Map options = [:]
    options.args          = args.args ?: ''
    options.args2         = args.args2 ?: ''
    options.publish_by_id = args.publish_by_id ?: false
    options.publish_dir   = args.publish_dir ?: ''
    options.publish_files = args.publish_files
    options.suffix        = args.suffix ?: ''
    return options
}

/*
 * Tidy up and join elements of a list to return a path string
 */
def getPathFromList(path_list) {
    def paths = path_list.findAll { item -> !item?.trim().isEmpty() }  // Remove empty entries
    paths = paths.collect { it.trim().replaceAll("^[/]+|[/]+\$", "") } // Trim whitespace and trailing slashes
    return paths.join('/')
}

/*
 * Function to save/publish module results
 */
def saveFiles(Map args) {
    if (!args.filename.endsWith('.version.txt')) {
        def ioptions = initOptions(args.options)
        def path_list = [ ioptions.publish_dir ?: args.publish_dir ]
        if (ioptions.publish_by_id) {
            path_list.add(args.publish_id)
        }
        if (ioptions.publish_files instanceof Map) {
            for (ext in ioptions.publish_files) {
                if (args.filename.endsWith(ext.key)) {
                    def ext_list = path_list.collect()
                    ext_list.add(ext.value)
                    return "${getPathFromList(ext_list)}/$args.filename"
                }
            }
        } else if (ioptions.publish_files == null) {
            return "${getPathFromList(path_list)}/$args.filename"
        }
    }
}
