#!/usr/bin/env nextflow

/*
================================================================================
                                  nf-core/sarek
================================================================================
Started March 2016.
Ported to nf-core May 2019.
--------------------------------------------------------------------------------
nf-core/sarek:
  An open-source analysis pipeline to detect germline or somatic variants
  from whole genome or targeted sequencing
--------------------------------------------------------------------------------
 @Homepage
 https://nf-co.re/sarek
--------------------------------------------------------------------------------
 @Documentation
 https://nf-co.re/sarek/docs
--------------------------------------------------------------------------------
*/

nextflow.preview.dsl = 2

/*
 * Print help message if required
 */
if (params.help) {
    def command = "nextflow run nf-core/sarek --input sample.tsv -profile docker"
    log.info Headers.nf_core(workflow, params.monochrome_logs)
    log.info Schema.params_help("$baseDir/nextflow_schema.json", command)
    exit 0
}

/*
================================================================================
                        INCLUDE SAREK FUNCTIONS
================================================================================
*/
include {
    check_parameter_existence;
    check_parameter_list;
    define_step_list;
    define_tool_list;
    extract_bam;
    extract_fastq;
    extract_fastq_from_dir;
    extract_infos;
    has_extension
} from './modules/local/functions'

/*
================================================================================
                         SET UP CONFIGURATION VARIABLES
================================================================================
*/

/*
 * Check parameters
 */

Checks.aws_batch(workflow, params)     // Check AWS batch settings
Checks.hostname(workflow, params, log) // Check the hostnames against configured profiles

/*
 * MultiQC
 * Stage config files
 */

multiqc_config = file("$baseDir/assets/multiqc_config.yaml", checkIfExists: true)
multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config, checkIfExists: true) : Channel.empty()
output_docs = file("$baseDir/docs/output.md", checkIfExists: true)
output_docs_images = file("$baseDir/docs/images/", checkIfExists: true)

// // Check if genome exists in the config file
// if (params.genomes && !params.genomes.containsKey(params.genome) && !params.igenomes_ignore) {
//     exit 1, "The provided genome '${params.genome}' is not available in the iGenomes file. Currently the available genomes are ${params.genomes.keySet().join(", ")}"
// } else if (params.genomes && !params.genomes.containsKey(params.genome) && params.igenomes_ignore) {
//     exit 1, "The provided genome '${params.genome}' is not available in the genomes file. Currently the available genomes are ${params.genomes.keySet().join(", ")}"
// }

step_list = define_step_list()
step = params.step ? params.step.toLowerCase().replaceAll('-', '').replaceAll('_', '') : ''

// Handle deprecation
if (step == 'preprocessing') step = 'mapping'

if (step.contains(',')) exit 1, 'You can choose only one step, see --help for more information'
if (!check_parameter_existence(step, step_list)) exit 1, "Unknown step ${step}, see --help for more information"

tool_list = define_tool_list()
tools = params.tools ? params.tools.split(',').collect{it.trim().toLowerCase().replaceAll('-', '').replaceAll('_', '')} : []
if (step == 'controlfreec') tools = ['controlfreec']
if (!check_parameter_list(tools, tool_list)) exit 1, 'Unknown tool(s), see --help for more information'

// skip__qc_list = define_skip_qc_list()
// skipQC = params.skip_qc ? params.skip_qc == 'all' ? skip__qc_list : params.skip_qc.split(',').collect{it.trim().toLowerCase().replaceAll('-', '').replaceAll('_', '')} : []
// if (!check_parameter_list(skipQC, skip__qc_list)) exit 1, 'Unknown QC tool(s), see --help for more information'

// anno_list = define_anno_list()
// annotate_tools = params.annotate_tools ? params.annotate_tools.split(',').collect{it.trim().toLowerCase().replaceAll('-', '')} : []
// if (!check_parameter_list(annotate_tools,anno_list)) exit 1, 'Unknown tool(s) to annotate, see --help for more information'

// // Check parameters
// if ((params.ascat_ploidy && !params.ascat_purity) || (!params.ascat_ploidy && params.ascat_purity)) exit 1, 'Please specify both --ascat_purity and --ascat_ploidy, or none of them'
// if (params.cf_window && params.cf_coeff) exit 1, 'Please specify either --cf_window OR --cf_coeff, but not both of them'


// Handle input
tsv_path = null
if (params.input && (has_extension(params.input, "tsv") || has_extension(params.input, "vcf") || has_extension(params.input, "vcf.gz"))) tsv_path = params.input
if (params.input && (has_extension(params.input, "vcf") || has_extension(params.input, "vcf.gz"))) step = "annotate"

// save_bam_mapped = params.skip_markduplicates ? true : params.save_bam_mapped ? true : false

// If no input file specified, trying to get TSV files corresponding to step in the TSV directory
// only for steps preparerecalibration, recalibrate, variantcalling and controlfreec
if (!params.input && params.sentieon) {
    switch (step) {
        case 'mapping': break
        case 'recalibrate': tsv_path = "${params.outdir}/Preprocessing/TSV/sentieon_deduped.tsv"; break
        case 'variantcalling': tsv_path = "${params.outdir}/Preprocessing/TSV/sentieon_recalibrated.tsv"; break
        case 'annotate': break
        default: exit 1, "Unknown step ${step}"
    }
} else if (!params.input && !params.sentieon && !params.skip_markduplicates) {
    switch (step) {
        case 'mapping': break
        case 'preparerecalibration': tsv_path = "${params.outdir}/Preprocessing/TSV/duplicates_marked_no_table.tsv"; break
        case 'recalibrate': tsv_path = "${params.outdir}/Preprocessing/TSV/duplicates_marked.tsv"; break
        case 'variantcalling': tsv_path = "${params.outdir}/Preprocessing/TSV/recalibrated.tsv"; break
        case 'controlfreec': tsv_path = "${params.outdir}/VariantCalling/TSV/control-freec_mpileup.tsv"; break
        case 'annotate': break
        default: exit 1, "Unknown step ${step}"
    }
} else if (!params.input && !params.sentieon && params.skip_markduplicates) {
    switch (step) {
        case 'mapping': break
        case 'preparerecalibration': tsv_path = "${params.outdir}/Preprocessing/TSV/mapped.tsv"; break
        case 'recalibrate': tsv_path = "${params.outdir}/Preprocessing/TSV/mapped_no_duplicates_marked.tsv"; break
        case 'variantcalling': tsv_path = "${params.outdir}/Preprocessing/TSV/recalibrated.tsv"; break
        case 'controlfreec': tsv_path = "${params.outdir}/VariantCalling/TSV/control-freec_mpileup.tsv"; break
        case 'annotate': break
        default: exit 1, "Unknown step ${step}"
    }
}

input_sample = Channel.empty()
if (tsv_path) {
    tsv_file = file(tsv_path)
    switch (step) {
        case 'mapping': input_sample = extract_fastq(tsv_file); break
        case 'preparerecalibration': input_sample = extract_bam(tsv_file); break
        case 'recalibrate': input_sample = extract_recal(tsv_file); break
        case 'variantcalling': input_sample = extract_bam(tsv_file); break
        case 'controlfreec': input_sample = extract_pileup(tsv_file); break
        case 'annotate': break
        default: exit 1, "Unknown step ${step}"
    }
} else if (params.input && !has_extension(params.input, "tsv")) {
    log.info "No TSV file"
    if (step != 'mapping') exit 1, 'No step other than "mapping" supports a directory as an input'
    log.info "Reading ${params.input} directory"
    log.warn "[nf-core/sarek] in ${params.input} directory, all fastqs are assuming to be from the same sample, which is assumed to be a germline one"
    input_sample = extract_fastq_from_dir(params.input)
    (input_sample, fastq_tmp) = input_sample.into(2)
    fastq_tmp.toList().subscribe onNext: {
        if (it.size() == 0) exit 1, "No FASTQ files found in --input directory '${params.input}'"
    }
    tsv_file = params.input  // used in the reports
} else if (tsv_path && step == 'annotate') {
    log.info "Annotating ${tsv_path}"
} else if (step == 'annotate') {
    log.info "Trying automatic annotation on files in the VariantCalling/ directory"
} else exit 1, 'No sample were defined, see --help'

(gender_map, status_map, input_sample) = extract_infos(input_sample)

// input_sample.dump(tag: 'input sample')

/*
================================================================================
                               CHECKING REFERENCES
================================================================================
*/

// Initialize each params in params.genomes, catch the command line first if it was defined
// params.fasta has to be the first one
params.fasta = params.genome && !('annotate' in step) ? params.genomes[params.genome].fasta ?: null : null

// The rest can be sorted
params.ac_loci = params.genome && 'ascat' in tools ? params.genomes[params.genome].ac_loci ?: null : null
params.ac_loci_gc = params.genome && 'ascat' in tools ? params.genomes[params.genome].ac_loci_gc ?: null : null
params.bwa = params.genome && params.fasta && 'mapping' in step ? params.genomes[params.genome].bwa ?: null : null
params.chr_dir = params.genome && 'controlfreec' in tools ? params.genomes[params.genome].chr_dir ?: null : null
params.chr_length = params.genome && 'controlfreec' in tools ? params.genomes[params.genome].chr_length ?: null : null
params.dbsnp = params.genome && ('mapping' in step || 'preparerecalibration' in step || 'controlfreec' in tools || 'haplotypecaller' in tools || 'mutect2' in tools || params.sentieon) ? params.genomes[params.genome].dbsnp ?: null : null
params.dbsnp_index = params.genome && params.dbsnp ? params.genomes[params.genome].dbsnp_index ?: null : null
params.dict = params.genome && params.fasta ? params.genomes[params.genome].dict ?: null : null
params.fasta_fai = params.genome && params.fasta ? params.genomes[params.genome].fasta_fai ?: null : null
params.germline_resource = params.genome && 'mutect2' in tools ? params.genomes[params.genome].germline_resource ?: null : null
params.germline_resource_index = params.genome && params.germline_resource ? params.genomes[params.genome].germline_resource_index ?: null : null
params.intervals = params.genome && !('annotate' in step) ? params.genomes[params.genome].intervals ?: null : null
params.known_indels = params.genome && ('mapping' in step || 'preparerecalibration' in step) ? params.genomes[params.genome].known_indels ?: null : null
params.known_indels_index = params.genome && params.known_indels ? params.genomes[params.genome].known_indels_index ?: null : null
params.mappability = params.genome && 'controlfreec' in tools ? params.genomes[params.genome].mappability ?: null : null
params.snpeff_db = params.genome && 'snpeff' in tools ? params.genomes[params.genome].snpeff_db ?: null : null
params.species = params.genome && 'vep' in tools ? params.genomes[params.genome].species ?: null : null
params.vep_cache_version = params.genome && 'vep' in tools ? params.genomes[params.genome].vep_cache_version ?: null : null

// Initialize channels based on params
fasta = params.fasta && !('annotate' in step) ? Channel.value(file(params.fasta)) : "null"
dbsnp = params.dbsnp && ('mapping' in step || 'preparerecalibration' in step || 'controlfreec' in tools || 'haplotypecaller' in tools || 'mutect2' in tools || params.sentieon) ? Channel.value(file(params.dbsnp)) : "null"
germline_resource = params.germline_resource && 'mutect2' in tools ? Channel.value(file(params.germline_resource)) : "null"
known_indels = params.known_indels && ('mapping' in step || 'preparerecalibration' in step) ? Channel.value(file(params.known_indels)) : "null"
pon = params.pon ? Channel.value(file(params.pon)) : "null"

ch_ac_loci = params.ac_loci && 'ascat' in tools ? Channel.value(file(params.ac_loci)) : "null"
ch_ac_loci_gc = params.ac_loci_gc && 'ascat' in tools ? Channel.value(file(params.ac_loci_gc)) : "null"
ch_chr_dir = params.chr_dir && 'controlfreec' in tools ? Channel.value(file(params.chr_dir)) : "null"
ch_chr_length = params.chr_length && 'controlfreec' in tools ? Channel.value(file(params.chr_length)) : "null"
fai = params.fasta_fai && !('annotate' in step) ? Channel.value(file(params.fasta_fai)) : "null"
intervals = params.intervals && !params.no_intervals && !('annotate' in step) ? Channel.value(file(params.intervals)) : "null"
ch_mappability = params.mappability && 'controlfreec' in tools ? Channel.value(file(params.mappability)) : "null"

ch_snpeff_cache = params.snpeff_cache ? Channel.value(file(params.snpeff_cache)) : "null"
ch_snpeff_db = params.snpeff_db ? Channel.value(params.snpeff_db) : "null"
ch_vep_cache_version = params.vep_cache_version ? Channel.value(params.vep_cache_version) : "null"
ch_vep_cache = params.vep_cache ? Channel.value(file(params.vep_cache)) : "null"

// Optional files, not defined within the params.genomes[params.genome] scope
ch_cadd_indels = params.cadd_indels ? Channel.value(file(params.cadd_indels)) : "null"
ch_cadd_indels_tbi = params.cadd_indels_tbi ? Channel.value(file(params.cadd_indels_tbi)) : "null"
ch_cadd_wg_snvs = params.cadd_wg_snvs ? Channel.value(file(params.cadd_wg_snvs)) : "null"
ch_cadd_wg_snvs_tbi = params.cadd_wg_snvs_tbi ? Channel.value(file(params.cadd_wg_snvs_tbi)) : "null"
ch_target_bed = params.target_bed ? Channel.value(file(params.target_bed)) : "null"

/*
================================================================================
                                PRINTING SUMMARY
================================================================================
*/

// Has the run name been specified by the user?
// This has the bonus effect of catching both -name and --name
run_name = params.name
if (!(workflow.runName ==~ /[a-z]+_[a-z]+/)) {
    run_name = workflow.runName
}
summary = Schema.params_summary(workflow, params, run_name)
log.info Headers.nf_core(workflow, params.monochrome_logs)
log.info summary.collect { k,v -> "${k.padRight(20)}: $v" }.join("\n")
log.info "-\033[2m----------------------------------------------------\033[0m-"

// params summary for MultiQC
workflow_summary = Schema.params_mqc_summary(summary)
workflow_summary = Channel.value(workflow_summary)

if ('mutect2' in tools && !(params.pon)) log.warn "[nf-core/sarek] Mutect2 was requested, but as no panel of normals were given, results will not be optimal"
if (params.sentieon) log.warn "[nf-core/sarek] Sentieon will be used, only works if Sentieon is available where nf-core/sarek is run"

/*
================================================================================
                        INCLUDE LOCAL PIPELINE MODULES
================================================================================
*/

include { BWAMEM2_MEM }           from './modules/local/bwamem2_mem.nf' addParams(params)
include { CREATE_INTERVALS_BED }  from './modules/local/create_intervals_bed' addParams(params)
include { GET_SOFTWARE_VERSIONS } from './modules/local/get_software_versions' params(params)
include { OUTPUT_DOCUMENTATION }  from './modules/local/output_documentation' params(params)
include { TRIM_GALORE }           from './modules/local/trim_galore.nf' addParams(params)

/*
================================================================================
                        INCLUDE nf-core PIPELINE MODULES
================================================================================
*/

include { FASTQC }  from './modules/nf-core/fastqc' params(params)
include { MULTIQC } from './modules/nf-core/multiqc' params(params)

// PREPARING CHANNELS FOR PREPROCESSING AND QC

// input_bam = Channel.empty()
// input_pair_reads = Channel.empty()

// if (step in ['preparerecalibration', 'recalibrate', 'variantcalling', 'controlfreec', 'annotate']) {
//     input_bam.close()
//     input_pair_reads.close()
// } else input_sample.branch(input_pair_reads, input_bam) {has_extension(it[3], "bam") ? 1 : 0}

// (input_bam, input_bam_fastqc) = input_bam.into(2)

// // Removing inputFile2 which is null in case of uBAM
// input_bam_fastqc = input_bam_fastqc.map {
//     idPatient, idSample, idRun, inputFile1, inputFile2 ->
//     [idPatient, idSample, idRun, inputFile1]
// }

// if (params.split_fastq){
//    input_pair_reads = input_pair_reads
//        // newly splitfastq are named based on split, so the name is easier to catch
//        .splitFastq(by: params.split_fastq, compress:true, file:"split", pe:true)
//        .map {idPatient, idSample, idRun, reads1, reads2 ->
//            // The split fastq read1 is the 4th element (indexed 3) its name is split_3
//            // The split fastq read2's name is split_4
//            // It's followed by which split it's acutally based on the mother fastq file
//            // Index start at 1
//            // Extracting the index to get a new IdRun
//            splitIndex = reads1.fileName.toString().minus("split_3.").minus(".gz")
//            newIdRun = idRun + "_" + splitIndex
//            // Giving the files a new nice name
//            newReads1 = file("${idSample}_${newIdRun}_R1.fastq.gz")
//            newReads2 = file("${idSample}_${newIdRun}_R2.fastq.gz")
//            [idPatient, idSample, newIdRun, reads1, reads2]}
//}

// input_pair_reads.dump(tag:'INPUT')

// (input_pair_reads, input_pair_readstrimgalore, input_pair_readsfastqc) = input_pair_reads.into(3)


/*
================================================================================
                        RUN THE WORKFLOW
================================================================================
*/

include { BUILD_INDICES } from './modules/subworkflows/build_indices' addParams(params)

fasta.dump(tag: 'fasta')

workflow {

    BUILD_INDICES(
        dbsnp,
        fasta,
        germline_resource,
        known_indels,
        pon,
        step)

    bwa = params.bwa ? Channel.value(file(params.bwa)) : BUILD_INDICES.out.bwa_built
    dict = params.dict ? Channel.value(file(params.dict)) : BUILD_INDICES.out.dictBuilt
    fai = params.fasta_fai ? Channel.value(file(params.fasta_fai)) : BUILD_INDICES.out.fai_built
    dbsnp_tbi = params.dbsnp ? params.dbsnp_index ? Channel.value(file(params.dbsnp_index)) : BUILD_INDICES.out.dbsnp_tbi : "null"
    germline_resource_tbi = params.germline_resource ? params.germline_resource_index ? Channel.value(file(params.germline_resource_index)) : BUILD_INDICES.out.germline_resource_tbi : "null"
    known_indels_tbi = params.known_indels ? params.known_indels_index ? Channel.value(file(params.known_indels_index)) : BUILD_INDICES.out.known_indels_tbi.collect() : "null"
    pon_tbi = params.pon ? params.pon_index ? Channel.value(file(params.pon_index)) : BUILD_INDICES.out.pon_tbi : "null"
    intervals = params.no_intervals ? "null" : params.intervals && !('annotate' in step) ? Channel.value(file(params.intervals)) : BUILD_INDICES.out.intervalBuilt
    intervals.dump(tag: 'intervals')

    // PREPROCESSING
    if((!params.no_intervals) && step != 'annotate')
        CREATE_INTERVALS_BED(intervals)

    // BED INTERVAL CHANNEL TRANSFORMING
    ch_bed_intervals = CREATE_INTERVALS_BED.out
        .flatten()
        .map { intervalFile ->
            def duration = 0.0
            for (line in intervalFile.readLines()) {
                final fields = line.split('\t')
                if (fields.size() >= 5) duration += fields[4].toFloat()
                else {
                    start = fields[1].toInteger()
                    end = fields[2].toInteger()
                    duration += (end - start) / params.nucleotides_per_second
                }
            }
            [ duration, intervalFile]
        }.toSortedList({ a, b -> b[0] <=> a[0] })
        .flatten().collate(2)
        .map{duration, intervalFile -> intervalFile}
    ch_bed_intervals.dump(tag:'bedintervals')

    if (params.no_intervals && step != 'annotate') {
        file("${params.outdir}/no_intervals.bed").text = "no_intervals\n"
        ch_bed_intervals = Channel.from(file("${params.outdir}/no_intervals.bed"))
    }

    //if(!('fastqc' in skipQC))
    FASTQC(input_sample)

    if(params.trim_fastq) {
        TRIM_GALORE(input_sample)
        BWAMEM2_MEM(TRIM_GALORE.out.trimmed_reads, bwa, fasta, fai)
    }
    else {
        BWAMEM2_MEM(input_sample, bwa, fasta, fai)
    }

    OUTPUT_DOCUMENTATION(
        output_docs,
        output_docs_images)

    GET_SOFTWARE_VERSIONS()

    MULTIQC(
        FASTQC.out.ifEmpty([]),
        multiqc_config,
        multiqc_custom_config.ifEmpty([]),
        GET_SOFTWARE_VERSIONS.out.yml,
        TRIM_GALORE.out.report.ifEmpty([]),
        workflow_summary)
}

/*
================================================================================
                        SEND COMPLETION EMAIL
================================================================================
 */
workflow.onComplete {
    def multiqc_report = []
    Completion.email(workflow, params, summary, run_name, baseDir, multiqc_report, log)
    Completion.summary(workflow, params, log)
}

// /*
// ================================================================================
//                                   PREPROCESSING
// ================================================================================
// */

//   (intBaseRecalibrator, intApplyBQSR, intHaplotypeCaller, intFreebayesSingle, intMpileup, bedIntervals) = bedIntervals.into(6)


// // STEP 0.5: QC ON READS

// // TODO: Use only one process for FastQC for FASTQ files and uBAM files
// // FASTQ and uBAM files are renamed based on the sample name


// process FastQCBAM {
//     label 'FastQC'
//     label 'cpus_2'

//     tag "${idPatient}-${idRun}"

//     publishDir "${params.outdir}/Reports/${idSample}/FastQC/${idSample}_${idRun}", mode: params.publish_dir_mode

//     input:
//         set idPatient, idSample, idRun, file("${idSample}_${idRun}.bam") from input_bam_fastqc

//     output:
//         file("*.{html,zip}") into fastQCBAMReport

//     when: !('fastqc' in skipQC)

//     script:
//     """
//     fastqc -t 2 -q ${idSample}_${idRun}.bam
//     """
// }

// fastQCReport = fastQCFQReport.mix(fastQCBAMReport)

// if (!params.trim_fastq) input_pair_readstrimgalore.close()

// // STEP 1: MAPPING READS TO REFERENCE GENOME WITH BWA MEM

// if (params.trim_fastq) input_pair_reads = outputPairReadsTrimGalore
// else input_pair_reads = input_pair_reads.mix(input_bam)

// input_pair_reads = input_pair_reads.dump(tag:'INPUT')

// (input_pair_reads, input_pair_reads_sentieon) = input_pair_reads.into(2)
// if (params.sentieon) input_pair_reads.close()
// else input_pair_reads_sentieon.close()

// bamMapped = bamMapped.dump(tag:'Mapped BAM')
// // Sort BAM whether they are standalone or should be merged

// singleBam = Channel.create()
// multipleBam = Channel.create()
// bamMapped.groupTuple(by:[0, 1])
//     .choice(singleBam, multipleBam) {it[2].size() > 1 ? 1 : 0}
// singleBam = singleBam.map {
//     idPatient, idSample, idRun, bam ->
//     [idPatient, idSample, bam]
// }
// singleBam = singleBam.dump(tag:'Single BAM')

// // STEP 1': MAPPING READS TO REFERENCE GENOME WITH SENTIEON BWA MEM

// process Sentieon_MapReads {
//     label 'cpus_max'
//     label 'memory_max'
//     label 'sentieon'

//     tag "${idPatient}-${idRun}"

//     input:
//         set idPatient, idSample, idRun, file(inputFile1), file(inputFile2) from input_pair_reads_sentieon
//         file(bwaIndex) from bwa
//         file(fasta) from fasta
//         file(fastaFai) from fai

//     output:
//         set idPatient, idSample, idRun, file("${idSample}_${idRun}.bam") into bam_sentieon_mapped

//     when: params.sentieon

//     script:
//     // -K is an hidden option, used to fix the number of reads processed by bwa mem
//     // Chunk size can affect bwa results, if not specified,
//     // the number of threads can change which can give not deterministic result.
//     // cf https://github.com/CCDG/Pipeline-Standardization/blob/master/PipelineStandard.md
//     // and https://github.com/gatk-workflows/gatk4-data-processing/blob/8ffa26ff4580df4ac3a5aa9e272a4ff6bab44ba2/processing-for-variant-discovery-gatk4.b37.wgs.inputs.json#L29
//     CN = params.sequencing_center ? "CN:${params.sequencing_center}\\t" : ""
//     readGroup = "@RG\\tID:${idRun}\\t${CN}PU:${idRun}\\tSM:${idSample}\\tLB:${idSample}\\tPL:illumina"
//     // adjust mismatch penalty for tumor samples
//     status = status_map[idPatient, idSample]
//     extra = status == 1 ? "-B 3" : ""
//     """
//     sentieon bwa mem -K 100000000 -R \"${readGroup}\" ${extra} -t ${task.cpus} -M ${fasta} \
//     ${inputFile1} ${inputFile2} | \
//     sentieon util sort -r ${fasta} -o ${idSample}_${idRun}.bam -t ${task.cpus} --sam2bam -i -
//     """
// }

// bam_sentieon_mapped = bam_sentieon_mapped.dump(tag:'Sentieon Mapped BAM')
// // Sort BAM whether they are standalone or should be merged

// singleBamSentieon = Channel.create()
// multipleBamSentieon = Channel.create()
// bam_sentieon_mapped.groupTuple(by:[0, 1])
//     .choice(singleBamSentieon, multipleBamSentieon) {it[2].size() > 1 ? 1 : 0}
// singleBamSentieon = singleBamSentieon.map {
//     idPatient, idSample, idRun, bam ->
//     [idPatient, idSample, bam]
// }
// singleBamSentieon = singleBamSentieon.dump(tag:'Single BAM')

// // STEP 1.5: MERGING BAM FROM MULTIPLE LANES

// multipleBam = multipleBam.mix(multipleBamSentieon)

// process MergeBamMapped {
//     label 'cpus_8'

//     tag "${idPatient}-${idSample}"

//     input:
//         set idPatient, idSample, idRun, file(bam) from multipleBam

//     output:
//         set idPatient, idSample, file("${idSample}.bam") into bam_mapped_merged

//     script:
//     """
//     samtools merge --threads ${task.cpus} ${idSample}.bam ${bam}
//     """
// }

// bam_mapped_merged = bam_mapped_merged.dump(tag:'Merged BAM')

// bam_mapped_merged = bam_mapped_merged.mix(singleBam,singleBamSentieon)

// (bam_mapped_merged, bam_sentieon_mapped_merged) = bam_mapped_merged.into(2)

// if (!params.sentieon) bam_sentieon_mapped_merged.close()
// else bam_mapped_merged.close()

// bam_mapped_merged = bam_mapped_merged.dump(tag:'BAMs for MD')
// bam_sentieon_mapped_merged = bam_sentieon_mapped_merged.dump(tag:'Sentieon BAMs to Index')

// process IndexBamMergedForSentieon {
//     label 'cpus_8'

//     tag "${idPatient}-${idSample}"

//     input:
//         set idPatient, idSample, file("${idSample}.bam") from bam_sentieon_mapped_merged

//     output:
//         set idPatient, idSample, file("${idSample}.bam"), file("${idSample}.bam.bai") into bam_sentieon_mapped_merged_indexed

//     script:
//     """
//     samtools index ${idSample}.bam
//     """
// }

// (bam_mapped_merged, bam_mapped_merged_to_index) = bam_mapped_merged.into(2)

// process IndexBamFile {
//     label 'cpus_8'

//     tag "${idPatient}-${idSample}"

//     publishDir params.outdir, mode: params.publish_dir_mode,
//         saveAs: {
//             if (save_bam_mapped) "Preprocessing/${idSample}/Mapped/${it}"
//             else null
//         }

//     input:
//         set idPatient, idSample, file("${idSample}.bam") from bam_mapped_merged_to_index

//     output:
//         set idPatient, idSample, file("${idSample}.bam"), file("${idSample}.bam.bai") into bam_mapped_merged_indexed
//         set idPatient, idSample into tsv_bam_indexed

//     when: save_bam_mapped || !(params.known_indels)

//     script:
//     """
//     samtools index ${idSample}.bam
//     """
// }

// if (!save_bam_mapped) tsv_bam_indexed.close()

// (tsv_bam_indexed, tsv_bam_indexed_sample) = tsv_bam_indexed.into(2)

// // Creating a TSV file to restart from this step
// tsv_bam_indexed.map { idPatient, idSample ->
//     gender = gender_map[idPatient]
//     status = status_map[idPatient, idSample]
//     bam = "${params.outdir}/Preprocessing/${idSample}/Mapped/${idSample}.bam"
//     bai = "${params.outdir}/Preprocessing/${idSample}/Mapped/${idSample}.bam.bai"
//     "${idPatient}\t${gender}\t${status}\t${idSample}\t${bam}\t${bai}\n"
// }.collectFile(
//     name: 'mapped.tsv', sort: true, storeDir: "${params.outdir}/Preprocessing/TSV"
// )

// tsv_bam_indexed_sample
//     .collectFile(storeDir: "${params.outdir}/Preprocessing/TSV") { idPatient, idSample ->
//         status = status_map[idPatient, idSample]
//         gender = gender_map[idPatient]
//         bam = "${params.outdir}/Preprocessing/${idSample}/Mapped/${idSample}.bam"
//         bai = "${params.outdir}/Preprocessing/${idSample}/Mapped/${idSample}.bam.bai"
//         ["mapped_${idSample}.tsv", "${idPatient}\t${gender}\t${status}\t${idSample}\t${bam}\t${bai}\n"]
// }
// // STEP 2: MARKING DUPLICATES

// process MarkDuplicates {
//     label 'cpus_16'

//     tag "${idPatient}-${idSample}"

//     publishDir params.outdir, mode: params.publish_dir_mode,
//         saveAs: {
//             if (it == "${idSample}.bam.metrics") "Reports/${idSample}/MarkDuplicates/${it}"
//             else "Preprocessing/${idSample}/DuplicatesMarked/${it}"
//         }

//     input:
//         set idPatient, idSample, file("${idSample}.bam") from bam_mapped_merged

//     output:
//         set idPatient, idSample, file("${idSample}.md.bam"), file("${idSample}.md.bam.bai") into bam_duplicates_marked
//         set idPatient, idSample into tsv_bam_duplicates_marked
//         file ("${idSample}.bam.metrics") optional true into duplicates_marked_report

//     when: !(params.skip_markduplicates)

//     script:
//     markdup_java_options = task.memory.toGiga() > 8 ? params.markdup_java_options : "\"-Xms" +  (task.memory.toGiga() / 2).trunc() + "g -Xmx" + (task.memory.toGiga() - 1) + "g\""
//     metrics = 'markduplicates' in skipQC ? '' : "-M ${idSample}.bam.metrics"
//     if (params.no_gatk_spark)
//     """
//     gatk --java-options ${markdup_java_options} \
//         MarkDuplicates \
//         --MAX_RECORDS_IN_RAM 50000 \
//         --INPUT ${idSample}.bam \
//         --METRICS_FILE ${idSample}.bam.metrics \
//         --TMP_DIR . \
//         --ASSUME_SORT_ORDER coordinate \
//         --CREATE_INDEX true \
//         --OUTPUT ${idSample}.md.bam

//     mv ${idSample}.md.bai ${idSample}.md.bam.bai
//     """
//     else
//     """
//     gatk --java-options ${markdup_java_options} \
//         MarkDuplicatesSpark \
//         -I ${idSample}.bam \
//         -O ${idSample}.md.bam \
//         ${metrics} \
//         --tmp-dir . \
//         --create-output-bam-index true \
//         --spark-master local[${task.cpus}]
//     """
// }

// (tsv_bam_duplicates_marked, tsv_bam_duplicates_marked_sample) = tsv_bam_duplicates_marked.into(2)

// // Creating a TSV file to restart from this step
// tsv_bam_duplicates_marked.map { idPatient, idSample ->
//     gender = gender_map[idPatient]
//     status = status_map[idPatient, idSample]
//     bam = "${params.outdir}/Preprocessing/${idSample}/DuplicatesMarked/${idSample}.md.bam"
//     bai = "${params.outdir}/Preprocessing/${idSample}/DuplicatesMarked/${idSample}.md.bam.bai"
//     "${idPatient}\t${gender}\t${status}\t${idSample}\t${bam}\t${bai}\n"
// }.collectFile(
//     name: 'duplicates_marked_no_table.tsv', sort: true, storeDir: "${params.outdir}/Preprocessing/TSV"
// )

// tsv_bam_duplicates_marked_sample
//     .collectFile(storeDir: "${params.outdir}/Preprocessing/TSV") { idPatient, idSample ->
//         status = status_map[idPatient, idSample]
//         gender = gender_map[idPatient]
//         bam = "${params.outdir}/Preprocessing/${idSample}/DuplicatesMarked/${idSample}.md.bam"
//         bai = "${params.outdir}/Preprocessing/${idSample}/DuplicatesMarked/${idSample}.md.bam.bai"
//         ["duplicates_marked_no_table_${idSample}.tsv", "${idPatient}\t${gender}\t${status}\t${idSample}\t${bam}\t${bai}\n"]
// }

// if ('markduplicates' in skipQC) duplicates_marked_report.close()

// if (step == 'preparerecalibration') bam_duplicates_marked = input_sample

// bam_duplicates_marked = bam_duplicates_marked.dump(tag:'MD BAM')
// duplicates_marked_report = duplicates_marked_report.dump(tag:'MD Report')

// if (params.skip_markduplicates) bam_duplicates_marked = bam_mapped_merged_indexed

// (bamMD, bamMDToJoin, bam_duplicates_marked) = bam_duplicates_marked.into(3)

// bamBaseRecalibrator = bamMD.combine(intBaseRecalibrator)

// bamBaseRecalibrator = bamBaseRecalibrator.dump(tag:'BAM FOR BASERECALIBRATOR')

// // STEP 2': SENTIEON DEDUP

// process Sentieon_Dedup {
//     label 'cpus_max'
//     label 'memory_max'
//     label 'sentieon'

//     tag "${idPatient}-${idSample}"

//     publishDir params.outdir, mode: params.publish_dir_mode,
//         saveAs: {
//             if (it == "${idSample}_*.txt" && 'sentieon' in skipQC) null
//             else if (it == "${idSample}_*.txt") "Reports/${idSample}/Sentieon/${it}"
//             else "Preprocessing/${idSample}/DedupedSentieon/${it}"
//         }

//     input:
//         set idPatient, idSample, file(bam), file(bai) from bam_sentieon_mapped_merged_indexed
//         file(fasta) from fasta
//         file(fastaFai) from fai

//     output:
//         set idPatient, idSample, file("${idSample}.deduped.bam"), file("${idSample}.deduped.bam.bai") into bam_sentieon_dedup

//     when: params.sentieon

//     script:
//     """
//     sentieon driver \
//         -t ${task.cpus} \
//         -i ${bam} \
//         -r ${fasta} \
//         --algo GCBias --summary ${idSample}_gc_summary.txt ${idSample}_gc_metric.txt \
//         --algo MeanQualityByCycle ${idSample}_mq_metric.txt \
//         --algo QualDistribution ${idSample}_qd_metric.txt \
//         --algo InsertSizeMetricAlgo ${idSample}_is_metric.txt  \
//         --algo AlignmentStat ${idSample}_aln_metric.txt

//     sentieon driver \
//         -t ${task.cpus} \
//         -i ${bam} \
//         --algo LocusCollector \
//         --fun score_info ${idSample}_score.gz

//     sentieon driver \
//         -t ${task.cpus} \
//         -i ${bam} \
//         --algo Dedup \
//         --rmdup \
//         --score_info ${idSample}_score.gz  \
//         --metrics ${idSample}_dedup_metric.txt ${idSample}.deduped.bam
//     """
// }

// // STEP 3: CREATING RECALIBRATION TABLES

// process BaseRecalibrator {
//     label 'cpus_1'

//     tag "${idPatient}-${idSample}-${intervalBed.baseName}"

//     input:
//         set idPatient, idSample, file(bam), file(bai), file(intervalBed) from bamBaseRecalibrator
//         file(dbsnp) from dbsnp
//         file(dbsnpIndex) from dbsnp_tbi
//         file(fasta) from fasta
//         file(dict) from dict
//         file(fastaFai) from fai
//         file(knownIndels) from known_indels
//         file(knownIndelsIndex) from known_indels_tbi

//     output:
//         set idPatient, idSample, file("${prefix}${idSample}.recal.table") into tableGatherBQSRReports
//         set idPatient, idSample into recalTableTSVnoInt

//     when: params.known_indels

//     script:
//     dbsnpOptions = params.dbsnp ? "--known-sites ${dbsnp}" : ""
//     knownOptions = params.known_indels ? knownIndels.collect{"--known-sites ${it}"}.join(' ') : ""
//     prefix = params.no_intervals ? "" : "${intervalBed.baseName}_"
//     intervalsOptions = params.no_intervals ? "" : "-L ${intervalBed}"
//     // TODO: --use-original-qualities ???
//     """
//     gatk --java-options -Xmx${task.memory.toGiga()}g \
//         BaseRecalibrator \
//         -I ${bam} \
//         -O ${prefix}${idSample}.recal.table \
//         --tmp-dir . \
//         -R ${fasta} \
//         ${intervalsOptions} \
//         ${dbsnpOptions} \
//         ${knownOptions} \
//         --verbosity INFO
//     """
// }

// if (!params.no_intervals) tableGatherBQSRReports = tableGatherBQSRReports.groupTuple(by:[0, 1])

// tableGatherBQSRReports = tableGatherBQSRReports.dump(tag:'BQSR REPORTS')

// if (params.no_intervals) {
//     (tableGatherBQSRReports, tableGatherBQSRReportsNoInt) = tableGatherBQSRReports.into(2)
//     recalTable = tableGatherBQSRReportsNoInt
// } else recalTableTSVnoInt.close()

// // STEP 3.5: MERGING RECALIBRATION TABLES

// process GatherBQSRReports {
//     label 'memory_singleCPU_2_task'
//     label 'cpus_2'

//     tag "${idPatient}-${idSample}"

//     publishDir params.outdir, mode: params.publish_dir_mode,
//         saveAs: {
//             if (it == "${idSample}.recal.table" && !params.skip_markduplicates) "Preprocessing/${idSample}/DuplicatesMarked/${it}"
//             else "Preprocessing/${idSample}/Mapped/${it}"
//         }

//     input:
//         set idPatient, idSample, file(recal) from tableGatherBQSRReports

//     output:
//         set idPatient, idSample, file("${idSample}.recal.table") into recalTable
//         file("${idSample}.recal.table") into baseRecalibratorReport
//         set idPatient, idSample into recalTableTSV

//     when: !(params.no_intervals)

//     script:
//     input = recal.collect{"-I ${it}"}.join(' ')
//     """
//     gatk --java-options -Xmx${task.memory.toGiga()}g \
//         GatherBQSRReports \
//         ${input} \
//         -O ${idSample}.recal.table \
//     """
// }

// if ('baserecalibrator' in skipQC) baseRecalibratorReport.close()

// recalTable = recalTable.dump(tag:'RECAL TABLE')

// (recalTableTSV, recalTableSampleTSV) = recalTableTSV.mix(recalTableTSVnoInt).into(2)

// // Create TSV files to restart from this step
// if (params.skip_markduplicates) {
//     recalTableTSV.map { idPatient, idSample ->
//         status = status_map[idPatient, idSample]
//         gender = gender_map[idPatient]
//         bam = "${params.outdir}/Preprocessing/${idSample}/Mapped/${idSample}.bam"
//         bai = "${params.outdir}/Preprocessing/${idSample}/Mapped/${idSample}.bam.bai"
//         recalTable = "${params.outdir}/Preprocessing/${idSample}/Mapped/${idSample}.recal.table"
//         "${idPatient}\t${gender}\t${status}\t${idSample}\t${bam}\t${bai}\t${recalTable}\n"
//     }.collectFile(
//         name: 'mapped_no_duplicates_marked.tsv', sort: true, storeDir: "${params.outdir}/Preprocessing/TSV"
//     )

//     recalTableSampleTSV
//         .collectFile(storeDir: "${params.outdir}/Preprocessing/TSV/") {
//             idPatient, idSample ->
//             status = status_map[idPatient, idSample]
//             gender = gender_map[idPatient]
//             bam = "${params.outdir}/Preprocessing/${idSample}/Mapped/${idSample}.bam"
//             bai = "${params.outdir}/Preprocessing/${idSample}/Mapped/${idSample}.bam.bai"
//             recalTable = "${params.outdir}/Preprocessing/${idSample}/Mapped/${idSample}.recal.table"
//             ["mapped_no_duplicates_marked_${idSample}.tsv", "${idPatient}\t${gender}\t${status}\t${idSample}\t${bam}\t${bai}\t${recalTable}\n"]
//     }
// } else {
//     recalTableTSV.map { idPatient, idSample ->
//     status = status_map[idPatient, idSample]
//     gender = gender_map[idPatient]
//     bam = "${params.outdir}/Preprocessing/${idSample}/DuplicatesMarked/${idSample}.md.bam"
//     bai = "${params.outdir}/Preprocessing/${idSample}/DuplicatesMarked/${idSample}.md.bam.bai"
//     recalTable = "${params.outdir}/Preprocessing/${idSample}/DuplicatesMarked/${idSample}.recal.table"

//         "${idPatient}\t${gender}\t${status}\t${idSample}\t${bam}\t${bai}\t${recalTable}\n"
//     }.collectFile(
//         name: 'duplicates_marked.tsv', sort: true, storeDir: "${params.outdir}/Preprocessing/TSV"
//     )

//     recalTableSampleTSV
//         .collectFile(storeDir: "${params.outdir}/Preprocessing/TSV/") {
//             idPatient, idSample ->
//             status = status_map[idPatient, idSample]
//             gender = gender_map[idPatient]
//             bam = "${params.outdir}/Preprocessing/${idSample}/DuplicatesMarked/${idSample}.md.bam"
//             bai = "${params.outdir}/Preprocessing/${idSample}/DuplicatesMarked/${idSample}.md.bam.bai"
//             recalTable = "${params.outdir}/Preprocessing/${idSample}/DuplicatesMarked/${idSample}.recal.table"
//             ["duplicates_marked_${idSample}.tsv", "${idPatient}\t${gender}\t${status}\t${idSample}\t${bam}\t${bai}\t${recalTable}\n"]
//     }
// }

// bamApplyBQSR = bamMDToJoin.join(recalTable, by:[0,1])

// if (step == 'recalibrate') bamApplyBQSR = input_sample

// bamApplyBQSR = bamApplyBQSR.dump(tag:'BAM + BAI + RECAL TABLE')

// bamApplyBQSR = bamApplyBQSR.combine(intApplyBQSR)

// bamApplyBQSR = bamApplyBQSR.dump(tag:'BAM + BAI + RECAL TABLE + INT')

// // STEP 4: RECALIBRATING

// process ApplyBQSR {
//     label 'memory_singleCPU_2_task'
//     label 'cpus_2'

//     tag "${idPatient}-${idSample}-${intervalBed.baseName}"

//     input:
//         set idPatient, idSample, file(bam), file(bai), file(recalibrationReport), file(intervalBed) from bamApplyBQSR
//         file(dict) from dict
//         file(fasta) from fasta
//         file(fastaFai) from fai

//     output:
//         set idPatient, idSample, file("${prefix}${idSample}.recal.bam") into bam_recalibrated_to_merge

//     script:
//     prefix = params.no_intervals ? "" : "${intervalBed.baseName}_"
//     intervalsOptions = params.no_intervals ? "" : "-L ${intervalBed}"
//     """
//     gatk --java-options -Xmx${task.memory.toGiga()}g \
//         ApplyBQSR \
//         -R ${fasta} \
//         --input ${bam} \
//         --output ${prefix}${idSample}.recal.bam \
//         ${intervalsOptions} \
//         --bqsr-recal-file ${recalibrationReport}
//     """
// }

// (bam_recalibrated_to_merge, bam_recalibrated_to_index) = bam_recalibrated_to_merge.groupTuple(by:[0, 1]).into(2)

// // STEP 4': SENTIEON BQSR

// bam_sentieon_dedup = bam_sentieon_dedup.dump(tag:'deduped.bam')

// process Sentieon_BQSR {
//     label 'cpus_max'
//     label 'memory_max'
//     label 'sentieon'

//     tag "${idPatient}-${idSample}"

//     publishDir params.outdir, mode: params.publish_dir_mode,
//         saveAs: {
//             if (it == "${idSample}_recal_result.csv" && 'sentieon' in skipQC) "Reports/${idSample}/Sentieon/${it}"
//             else "Preprocessing/${idSample}/RecalSentieon/${it}"
//         }

//     input:
//         set idPatient, idSample, file(bam), file(bai) from bam_sentieon_dedup
//         file(dbsnp) from dbsnp
//         file(dbsnpIndex) from dbsnp_tbi
//         file(fasta) from fasta
//         file(dict) from dict
//         file(fastaFai) from fai
//         file(knownIndels) from known_indels
//         file(knownIndelsIndex) from known_indels_tbi

//     output:
//         set idPatient, idSample, file("${idSample}.recal.bam"), file("${idSample}.recal.bam.bai") into bam_sentieon_recal
//         set idPatient, idSample, file(bam), file(bai), file("${idSample}.recal.table") into bam_sentieon_deduped_table
//         set idPatient, idSample into tsv_sentieon

//     when: params.sentieon

//     script:
//     known = knownIndels.collect{"--known-sites ${it}"}.join(' ')
//     """
//     sentieon driver  \
//         -t ${task.cpus} \
//         -r ${fasta} \
//         -i ${idSample}.deduped.bam \
//         --algo QualCal \
//         -k ${dbsnp} \
//         ${idSample}.recal.table

//     sentieon driver \
//         -t ${task.cpus} \
//         -r ${fasta} \
//         -i ${idSample}.deduped.bam \
//         -q ${idSample}.recal.table \
//         --algo QualCal \
//         -k ${dbsnp} \
//         ${idSample}.table.post \
//         --algo ReadWriter ${idSample}.recal.bam

//     sentieon driver \
//         -t ${task.cpus} \
//         --algo QualCal \
//         --plot \
//         --before ${idSample}.recal.table \
//         --after ${idSample}.table.post \
//         ${idSample}_recal_result.csv
//     """
// }

// (tsv_sentieon_deduped, tsv_sentieon_deduped_sample, tsv_sentieon_recal, tsv_sentieon_recal_sample) = tsv_sentieon.into(4)

// // Creating a TSV file to restart from this step
// tsv_sentieon_deduped.map { idPatient, idSample ->
//     gender = gender_map[idPatient]
//     status = status_map[idPatient, idSample]
//     bam = "${params.outdir}/Preprocessing/${idSample}/DedupedSentieon/${idSample}.deduped.bam"
//     bai = "${params.outdir}/Preprocessing/${idSample}/DedupedSentieon/${idSample}.deduped.bam.bai"
//     table = "${params.outdir}/Preprocessing/${idSample}/RecalSentieon/${idSample}.recal.table"
//     "${idPatient}\t${gender}\t${status}\t${idSample}\t${bam}\t${bai}\t${table}\n"
// }.collectFile(
//     name: 'sentieon_deduped.tsv', sort: true, storeDir: "${params.outdir}/Preprocessing/TSV"
// )

// tsv_sentieon_deduped_sample
//     .collectFile(storeDir: "${params.outdir}/Preprocessing/TSV") { idPatient, idSample ->
//         status = status_map[idPatient, idSample]
//         gender = gender_map[idPatient]
//         bam = "${params.outdir}/Preprocessing/${idSample}/DedupedSentieon/${idSample}.deduped.bam"
//         bai = "${params.outdir}/Preprocessing/${idSample}/DedupedSentieon/${idSample}.deduped.bam.bai"
//         table = "${params.outdir}/Preprocessing/${idSample}/RecalSentieon/${idSample}.recal.table"
//         ["sentieon_deduped_${idSample}.tsv", "${idPatient}\t${gender}\t${status}\t${idSample}\t${bam}\t${bai}\t${table}\n"]
// }

// // Creating a TSV file to restart from this step
// tsv_sentieon_recal.map { idPatient, idSample ->
//     gender = gender_map[idPatient]
//     status = status_map[idPatient, idSample]
//     bam = "${params.outdir}/Preprocessing/${idSample}/RecalSentieon/${idSample}.recal.bam"
//     bai = "${params.outdir}/Preprocessing/${idSample}/RecalSentieon/${idSample}.recal.bam.bai"
//     "${idPatient}\t${gender}\t${status}\t${idSample}\t${bam}\t${bai}\n"
// }.collectFile(
//     name: 'sentieon_recalibrated.tsv', sort: true, storeDir: "${params.outdir}/Preprocessing/TSV"
// )

// tsv_sentieon_recal_sample
//     .collectFile(storeDir: "${params.outdir}/Preprocessing/TSV") { idPatient, idSample ->
//         status = status_map[idPatient, idSample]
//         gender = gender_map[idPatient]
//         bam = "${params.outdir}/Preprocessing/${idSample}/RecalSentieon/${idSample}.recal.bam"
//         bai = "${params.outdir}/Preprocessing/${idSample}/RecalSentieon/${idSample}.recal.bam.bai"
//         ["sentieon_recalibrated_${idSample}.tsv", "${idPatient}\t${gender}\t${status}\t${idSample}\t${bam}\t${bai}\n"]
// }

// // STEP 4.5: MERGING THE RECALIBRATED BAM FILES

// process MergeBamRecal {
//     label 'cpus_8'

//     tag "${idPatient}-${idSample}"

//     publishDir "${params.outdir}/Preprocessing/${idSample}/Recalibrated", mode: params.publish_dir_mode

//     input:
//         set idPatient, idSample, file(bam) from bam_recalibrated_to_merge

//     output:
//         set idPatient, idSample, file("${idSample}.recal.bam"), file("${idSample}.recal.bam.bai") into bam_recalibrated
//         set idPatient, idSample, file("${idSample}.recal.bam") into bam_recalibrated_qc
//         set idPatient, idSample into tsv_bam_recalibrated

//     when: !(params.no_intervals)

//     script:
//     """
//     samtools merge --threads ${task.cpus} ${idSample}.recal.bam ${bam}
//     samtools index ${idSample}.recal.bam
//     """
// }

// // STEP 4.5': INDEXING THE RECALIBRATED BAM FILES

// process IndexBamRecal {
//     label 'cpus_8'

//     tag "${idPatient}-${idSample}"

//     publishDir "${params.outdir}/Preprocessing/${idSample}/Recalibrated", mode: params.publish_dir_mode

//     input:
//         set idPatient, idSample, file("${idSample}.recal.bam") from bam_recalibrated_to_index

//     output:
//         set idPatient, idSample, file("${idSample}.recal.bam"), file("${idSample}.recal.bam.bai") into bam_recalibrated_indexed
//         set idPatient, idSample, file("${idSample}.recal.bam") into bam_recalibrated_no_int_qc
//         set idPatient, idSample into tsv_bam_recalibrated_no_int

//     when: params.no_intervals

//     script:
//     """
//     samtools index ${idSample}.recal.bam
//     """
// }

// bam_recalibrated = bam_recalibrated.mix(bam_recalibrated_indexed)
// bam_recalibrated_qc = bam_recalibrated_qc.mix(bam_recalibrated_no_int_qc)
// tsv_bam_recalibrated = tsv_bam_recalibrated.mix(tsv_bam_recalibrated_no_int)

// (bam_recalibrated_bamqc, bam_recalibrated_samtools_stats) = bam_recalibrated_qc.into(2)
// (tsv_bam_recalibrated, tsv_bam_recalibrated_sample) = tsv_bam_recalibrated.into(2)

// // Creating a TSV file to restart from this step
// tsv_bam_recalibrated.map { idPatient, idSample ->
//     gender = gender_map[idPatient]
//     status = status_map[idPatient, idSample]
//     bam = "${params.outdir}/Preprocessing/${idSample}/Recalibrated/${idSample}.recal.bam"
//     bai = "${params.outdir}/Preprocessing/${idSample}/Recalibrated/${idSample}.recal.bam.bai"
//     "${idPatient}\t${gender}\t${status}\t${idSample}\t${bam}\t${bai}\n"
// }.collectFile(
//     name: 'recalibrated.tsv', sort: true, storeDir: "${params.outdir}/Preprocessing/TSV"
// )

// tsv_bam_recalibrated_sample
//     .collectFile(storeDir: "${params.outdir}/Preprocessing/TSV") {
//         idPatient, idSample ->
//         status = status_map[idPatient, idSample]
//         gender = gender_map[idPatient]
//         bam = "${params.outdir}/Preprocessing/${idSample}/Recalibrated/${idSample}.recal.bam"
//         bai = "${params.outdir}/Preprocessing/${idSample}/Recalibrated/${idSample}.recal.bam.bai"
//         ["recalibrated_${idSample}.tsv", "${idPatient}\t${gender}\t${status}\t${idSample}\t${bam}\t${bai}\n"]
// }

// // STEP 5: QC

// process SamtoolsStats {
//     label 'cpus_2'

//     tag "${idPatient}-${idSample}"

//     publishDir "${params.outdir}/Reports/${idSample}/SamToolsStats", mode: params.publish_dir_mode

//     input:
//         set idPatient, idSample, file(bam) from bam_recalibrated_samtools_stats

//     output:
//         file ("${bam}.samtools.stats.out") into samtoolsStatsReport

//     when: !('samtools' in skipQC)

//     script:
//     """
//     samtools stats ${bam} > ${bam}.samtools.stats.out
//     """
// }

// samtoolsStatsReport = samtoolsStatsReport.dump(tag:'SAMTools')

// bamBamQC = bamMappedBamQC.mix(bam_recalibrated_bamqc)

// process BamQC {
//     label 'memory_max'
//     label 'cpus_16'

//     tag "${idPatient}-${idSample}"

//     publishDir "${params.outdir}/Reports/${idSample}/bamQC", mode: params.publish_dir_mode

//     input:
//         set idPatient, idSample, file(bam) from bamBamQC
//         file(targetBED) from ch_target_bed

//     output:
//         file("${bam.baseName}") into bamQCReport

//     when: !('bamqc' in skipQC)

//     script:
//     use_bed = params.target_bed ? "-gff ${targetBED}" : ''
//     """
//     qualimap --java-mem-size=${task.memory.toGiga()}G \
//         bamqc \
//         -bam ${bam} \
//         --paint-chromosome-limits \
//         --genome-gc-distr HUMAN \
//         $use_bed \
//         -nt ${task.cpus} \
//         -skip-duplicated \
//         --skip-dup-mode 0 \
//         -outdir ${bam.baseName} \
//         -outformat HTML
//     """
// }

// bamQCReport = bamQCReport.dump(tag:'BamQC')

// /*
// ================================================================================
//                             GERMLINE VARIANT CALLING
// ================================================================================
// */

// // When using sentieon for mapping, Channel bam_recalibrated is bam_sentieon_recal
// if (params.sentieon && step == 'mapping') bam_recalibrated = bam_sentieon_recal

// // When no knownIndels for mapping, Channel bam_recalibrated is bam_duplicates_marked
// if (!params.known_indels && step == 'mapping') bam_recalibrated = bam_duplicates_marked

// // When starting with variant calling, Channel bam_recalibrated is input_sample
// if (step == 'variantcalling') bam_recalibrated = input_sample

// bam_recalibrated = bam_recalibrated.dump(tag:'BAM for Variant Calling')

// // Here we have a recalibrated bam set
// // The TSV file is formatted like: "idPatient status idSample bamFile baiFile"
// // Manta will be run in Germline mode, or in Tumor mode depending on status
// // HaplotypeCaller, TIDDIT and Strelka will be run for Normal and Tumor samples

// (bamMantaSingle, bamStrelkaSingle, bamTIDDIT, bamFreebayesSingleNoIntervals, bamHaplotypeCallerNoIntervals, bamRecalAll) = bam_recalibrated.into(6)

// (bam_sentieon_DNAseq, bam_sentieon_DNAscope, bam_sentieon_all) = bam_sentieon_deduped_table.into(3)

// // To speed Variant Callers up we are chopping the reference into smaller pieces
// // Do variant calling by this intervals, and re-merge the VCFs

// bamHaplotypeCaller = bamHaplotypeCallerNoIntervals.spread(intHaplotypeCaller)
// bamFreebayesSingle = bamFreebayesSingleNoIntervals.spread(intFreebayesSingle)

// // STEP GATK HAPLOTYPECALLER.1

// process HaplotypeCaller {
//     label 'memory_singleCPU_task_sq'
//     label 'cpus_2'

//     tag "${idSample}-${intervalBed.baseName}"

//     input:
//         set idPatient, idSample, file(bam), file(bai), file(intervalBed) from bamHaplotypeCaller
//         file(dbsnp) from dbsnp
//         file(dbsnpIndex) from dbsnp_tbi
//         file(dict) from dict
//         file(fasta) from fasta
//         file(fastaFai) from fai

//     output:
//         set val("HaplotypeCallerGVCF"), idPatient, idSample, file("${intervalBed.baseName}_${idSample}.g.vcf") into gvcfHaplotypeCaller
//         set idPatient, idSample, file(intervalBed), file("${intervalBed.baseName}_${idSample}.g.vcf") into gvcfGenotypeGVCFs

//     when: 'haplotypecaller' in tools

//     script:
//     intervalsOptions = params.no_intervals ? "" : "-L ${intervalBed}"
//     dbsnpOptions = params.dbsnp ? "--D ${dbsnp}" : ""
//     """
//     gatk --java-options "-Xmx${task.memory.toGiga()}g -Xms6000m -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10" \
//         HaplotypeCaller \
//         -R ${fasta} \
//         -I ${bam} \
//         ${intervalsOptions} \
//         ${dbsnpOptions} \
//         -O ${intervalBed.baseName}_${idSample}.g.vcf \
//         -ERC GVCF
//     """
// }

// gvcfHaplotypeCaller = gvcfHaplotypeCaller.groupTuple(by:[0, 1, 2])

// if (params.no_gvcf) gvcfHaplotypeCaller.close()
// else gvcfHaplotypeCaller = gvcfHaplotypeCaller.dump(tag:'GVCF HaplotypeCaller')

// // STEP GATK HAPLOTYPECALLER.2

// process GenotypeGVCFs {
//     tag "${idSample}-${intervalBed.baseName}"

//     input:
//         set idPatient, idSample, file(intervalBed), file(gvcf) from gvcfGenotypeGVCFs
//         file(dbsnp) from dbsnp
//         file(dbsnpIndex) from dbsnp_tbi
//         file(dict) from dict
//         file(fasta) from fasta
//         file(fastaFai) from fai

//     output:
//     set val("HaplotypeCaller"), idPatient, idSample, file("${intervalBed.baseName}_${idSample}.vcf") into vcfGenotypeGVCFs

//     when: 'haplotypecaller' in tools

//     script:
//     // Using -L is important for speed and we have to index the interval files also
//     intervalsOptions = params.no_intervals ? "" : "-L ${intervalBed}"
//     dbsnpOptions = params.dbsnp ? "--D ${dbsnp}" : ""
//     """
//     gatk --java-options -Xmx${task.memory.toGiga()}g \
//         IndexFeatureFile \
//         -I ${gvcf}

//     gatk --java-options -Xmx${task.memory.toGiga()}g \
//         GenotypeGVCFs \
//         -R ${fasta} \
//         ${intervalsOptions} \
//         ${dbsnpOptions} \
//         -V ${gvcf} \
//         -O ${intervalBed.baseName}_${idSample}.vcf
//     """
// }

// vcfGenotypeGVCFs = vcfGenotypeGVCFs.groupTuple(by:[0, 1, 2])

// // STEP SENTIEON DNAseq

// process Sentieon_DNAseq {
//     label 'cpus_max'
//     label 'memory_max'
//     label 'sentieon'

//     tag "${idSample}"

//     input:
//         set idPatient, idSample, file(bam), file(bai), file(recal) from bam_sentieon_DNAseq
//         file(dbsnp) from dbsnp
//         file(dbsnpIndex) from dbsnp_tbi
//         file(fasta) from fasta
//         file(fastaFai) from fai

//     output:
//     set val("SentieonDNAseq"), idPatient, idSample, file("DNAseq_${idSample}.vcf") into vcf_sentieon_DNAseq

//     when: 'dnaseq' in tools && params.sentieon

//     script:
//     """
//     sentieon driver \
//         -t ${task.cpus} \
//         -r ${fasta} \
//         -i ${bam} \
//         -q ${recal} \
//         --algo Haplotyper \
//         -d ${dbsnp} \
//         DNAseq_${idSample}.vcf
//     """
// }

// vcf_sentieon_DNAseq = vcf_sentieon_DNAseq.dump(tag:'sentieon DNAseq')

// // STEP SENTIEON DNAscope

// process Sentieon_DNAscope {
//     label 'cpus_max'
//     label 'memory_max'
//     label 'sentieon'

//     tag "${idSample}"

//     input:
//         set idPatient, idSample, file(bam), file(bai), file(recal) from bam_sentieon_DNAscope
//         file(dbsnp) from dbsnp
//         file(dbsnpIndex) from dbsnp_tbi
//         file(fasta) from fasta
//         file(fastaFai) from fai

//     output:
//     set val("SentieonDNAscope"), idPatient, idSample, file("DNAscope_${idSample}.vcf") into vcf_sentieon_DNAscope
//     set val("SentieonDNAscope"), idPatient, idSample, file("DNAscope_SV_${idSample}.vcf") into vcf_sentieon_DNAscope_SV

//     when: 'dnascope' in tools && params.sentieon

//     script:
//     """
//     sentieon driver \
//         -t ${task.cpus} \
//         -r ${fasta} \
//         -i ${bam} \
//         -q ${recal} \
//         --algo DNAscope \
//         -d ${dbsnp} \
//         DNAscope_${idSample}.vcf

//     sentieon driver \
//         -t ${task.cpus} \
//         -r ${fasta}\
//         -i ${bam} \
//         -q ${recal} \
//         --algo DNAscope \
//         --var_type bnd \
//         -d ${dbsnp} \
//         DNAscope_${idSample}.temp.vcf

//     sentieon driver \
//         -t ${task.cpus} \
//         -r ${fasta}\
//         -q ${recal} \
//         --algo SVSolver \
//         -v DNAscope_${idSample}.temp.vcf \
//         DNAscope_SV_${idSample}.vcf
//     """
// }

// vcf_sentieon_DNAscope = vcf_sentieon_DNAscope.dump(tag:'sentieon DNAscope')
// vcf_sentieon_DNAscope_SV = vcf_sentieon_DNAscope_SV.dump(tag:'sentieon DNAscope SV')

// // STEP STRELKA.1 - SINGLE MODE

// process StrelkaSingle {
//     label 'cpus_max'
//     label 'memory_max'

//     tag "${idSample}"

//     publishDir "${params.outdir}/VariantCalling/${idSample}/Strelka", mode: params.publish_dir_mode

//     input:
//         set idPatient, idSample, file(bam), file(bai) from bamStrelkaSingle
//         file(fasta) from fasta
//         file(fastaFai) from fai
//         file(targetBED) from ch_target_bed

//     output:
//         set val("Strelka"), idPatient, idSample, file("*.vcf.gz"), file("*.vcf.gz.tbi") into vcfStrelkaSingle

//     when: 'strelka' in tools

//     script:
//     beforeScript = params.target_bed ? "bgzip --threads ${task.cpus} -c ${targetBED} > call_targets.bed.gz ; tabix call_targets.bed.gz" : ""
//     options = params.target_bed ? "--exome --callRegions call_targets.bed.gz" : ""
//     """
//     ${beforeScript}
//     configureStrelkaGermlineWorkflow.py \
//         --bam ${bam} \
//         --referenceFasta ${fasta} \
//         ${options} \
//         --runDir Strelka

//     python Strelka/runWorkflow.py -m local -j ${task.cpus}

//     mv Strelka/results/variants/genome.*.vcf.gz \
//         Strelka_${idSample}_genome.vcf.gz
//     mv Strelka/results/variants/genome.*.vcf.gz.tbi \
//         Strelka_${idSample}_genome.vcf.gz.tbi
//     mv Strelka/results/variants/variants.vcf.gz \
//         Strelka_${idSample}_variants.vcf.gz
//     mv Strelka/results/variants/variants.vcf.gz.tbi \
//         Strelka_${idSample}_variants.vcf.gz.tbi
//     """
// }

// vcfStrelkaSingle = vcfStrelkaSingle.dump(tag:'Strelka - Single Mode')

// // STEP MANTA.1 - SINGLE MODE

// process MantaSingle {
//     label 'cpus_max'
//     label 'memory_max'

//     tag "${idSample}"

//     publishDir "${params.outdir}/VariantCalling/${idSample}/Manta", mode: params.publish_dir_mode

//     input:
//         set idPatient, idSample, file(bam), file(bai) from bamMantaSingle
//         file(fasta) from fasta
//         file(fastaFai) from fai
//         file(targetBED) from ch_target_bed

//     output:
//         set val("Manta"), idPatient, idSample, file("*.vcf.gz"), file("*.vcf.gz.tbi") into vcfMantaSingle

//     when: 'manta' in tools

//     script:
//     beforeScript = params.target_bed ? "bgzip --threads ${task.cpus} -c ${targetBED} > call_targets.bed.gz ; tabix call_targets.bed.gz" : ""
//     options = params.target_bed ? "--exome --callRegions call_targets.bed.gz" : ""
//     status = status_map[idPatient, idSample]
//     input_bam = status == 0 ? "--bam" : "--tumorBam"
//     vcftype = status == 0 ? "diploid" : "tumor"
//     """
//     ${beforeScript}
//     configManta.py \
//         ${input_bam} ${bam} \
//         --reference ${fasta} \
//         ${options} \
//         --runDir Manta

//     python Manta/runWorkflow.py -m local -j ${task.cpus}

//     mv Manta/results/variants/candidateSmallIndels.vcf.gz \
//         Manta_${idSample}.candidateSmallIndels.vcf.gz
//     mv Manta/results/variants/candidateSmallIndels.vcf.gz.tbi \
//         Manta_${idSample}.candidateSmallIndels.vcf.gz.tbi
//     mv Manta/results/variants/candidateSV.vcf.gz \
//         Manta_${idSample}.candidateSV.vcf.gz
//     mv Manta/results/variants/candidateSV.vcf.gz.tbi \
//         Manta_${idSample}.candidateSV.vcf.gz.tbi
//     mv Manta/results/variants/${vcftype}SV.vcf.gz \
//         Manta_${idSample}.${vcftype}SV.vcf.gz
//     mv Manta/results/variants/${vcftype}SV.vcf.gz.tbi \
//         Manta_${idSample}.${vcftype}SV.vcf.gz.tbi
//     """
// }

// vcfMantaSingle = vcfMantaSingle.dump(tag:'Single Manta')

// // STEP TIDDIT

// process TIDDIT {
//     tag "${idSample}"

//     publishDir params.outdir, mode: params.publish_dir_mode,
//         saveAs: {
//             if (it == "TIDDIT_${idSample}.vcf") "VariantCalling/${idSample}/TIDDIT/${it}"
//             else "Reports/${idSample}/TIDDIT/${it}"
//         }

//     input:
//         set idPatient, idSample, file(bam), file(bai) from bamTIDDIT
//         file(fasta) from fasta
//         file(fastaFai) from fai

//     output:
//         set val("TIDDIT"), idPatient, idSample, file("*.vcf.gz"), file("*.tbi") into vcfTIDDIT
//         set file("TIDDIT_${idSample}.old.vcf"), file("TIDDIT_${idSample}.ploidy.tab"), file("TIDDIT_${idSample}.signals.tab"), file("TIDDIT_${idSample}.wig"), file("TIDDIT_${idSample}.gc.wig") into tidditOut

//     when: 'tiddit' in tools

//     script:
//     """
//     tiddit --sv -o TIDDIT_${idSample} --bam ${bam} --ref ${fasta}

//     mv TIDDIT_${idSample}.vcf TIDDIT_${idSample}.old.vcf

//     grep -E "#|PASS" TIDDIT_${idSample}.old.vcf > TIDDIT_${idSample}.vcf

//     bgzip --threads ${task.cpus} -c TIDDIT_${idSample}.vcf > TIDDIT_${idSample}.vcf.gz

//     tabix TIDDIT_${idSample}.vcf.gz
//     """
// }

// vcfTIDDIT = vcfTIDDIT.dump(tag:'TIDDIT')

// // STEP FREEBAYES SINGLE MODE

// process FreebayesSingle {
//     tag "${idSample}-${intervalBed.baseName}"

//     label 'cpus_1'

//     input:
//         set idPatient, idSample, file(bam), file(bai), file(intervalBed) from bamFreebayesSingle
//         file(fasta) from fasta
//         file(fastaFai) from ch_software_versions_yaml

//     output:
//         set val("FreeBayes"), idPatient, idSample, file("${intervalBed.baseName}_${idSample}.vcf") into vcfFreebayesSingle

//     when: 'freebayes' in tools

//     script:
//     intervalsOptions = params.no_intervals ? "" : "-t ${intervalBed}"
//     """
//     freebayes \
//         -f ${fasta} \
//         --min-alternate-fraction 0.1 \
//         --min-mapping-quality 1 \
//         ${intervalsOptions} \
//         ${bam} > ${intervalBed.baseName}_${idSample}.vcf
//     """
// }

// vcfFreebayesSingle = vcfFreebayesSingle.groupTuple(by: [0,1,2])

// /*
// ================================================================================
//                              SOMATIC VARIANT CALLING
// ================================================================================
// */
// // Ascat, pileup, pileups with no intervals, recalibrated BAMs
// (bamAscat, bamMpileup, bamMpileupNoInt, bamRecalAll) = bamRecalAll.into(4)

// // separate BAM by status
// bamNormal = Channel.create()
// bamTumor = Channel.create()

// bamRecalAll
//     .choice(bamTumor, bamNormal) {status_map[it[0], it[1]] == 0 ? 1 : 0}

// // Crossing Normal and Tumor to get a T/N pair for Somatic Variant Calling
// // Remapping channel to remove common key idPatient
// pairBam = bamNormal.cross(bamTumor).map {
//     normal, tumor ->
//     [normal[0], normal[1], normal[2], normal[3], tumor[1], tumor[2], tumor[3]]
// }

// pairBam = pairBam.dump(tag:'BAM Somatic Pair')

// // Manta, Strelka, Mutect2, MSIsensor
// (pairBamManta, pairBamStrelka, pairBamStrelkaBP, pairBamCalculateContamination, pairBamFilterMutect2, pairBamMsisensor, pairBamCNVkit, pairBam) = pairBam.into(8)

// // Making Pair Bam for Sention

// // separate BAM by status
// bam_sention_normal = Channel.create()
// bam_sentieon_tumor = Channel.create()

// bam_sentieon_all
//     .choice(bam_sentieon_tumor, bam_sention_normal) {status_map[it[0], it[1]] == 0 ? 1 : 0}

// // Crossing Normal and Tumor to get a T/N pair for Somatic Variant Calling
// // Remapping channel to remove common key idPatient

// bam_pair_sentieon_TNscope = bam_sention_normal.cross(bam_sentieon_tumor).map {
//     normal, tumor ->
//     [normal[0], normal[1], normal[2], normal[3], normal[4], tumor[1], tumor[2], tumor[3], tumor[4]]
// }

// intervalPairBam = pairBam.spread(bedIntervals)

// bamMpileup = bamMpileup.spread(intMpileup)

// // intervals for Mutect2 calls, FreeBayes and pileups for Mutect2 filtering
// (pairBamMutect2, pairBamFreeBayes, pairBamPileupSummaries) = intervalPairBam.into(3)

// // STEP FREEBAYES

// process FreeBayes {
//     tag "${idSampleTumor}_vs_${idSampleNormal}-${intervalBed.baseName}"

//     label 'cpus_1'

//     input:
//         set idPatient, idSampleNormal, file(bamNormal), file(baiNormal), idSampleTumor, file(bamTumor), file(baiTumor), file(intervalBed) from pairBamFreeBayes
//         file(fasta) from fasta
//         file(fastaFai) from fai

//     output:
//         set val("FreeBayes"), idPatient, val("${idSampleTumor}_vs_${idSampleNormal}"), file("${intervalBed.baseName}_${idSampleTumor}_vs_${idSampleNormal}.vcf") into vcfFreeBayes

//     when: 'freebayes' in tools

//     script:
//     intervalsOptions = params.no_intervals ? "" : "-t ${intervalBed}"
//     """
//     freebayes \
//         -f ${fasta} \
//         --pooled-continuous \
//         --pooled-discrete \
//         --genotype-qualities \
//         --report-genotype-likelihood-max \
//         --allele-balance-priors-off \
//         --min-alternate-fraction 0.03 \
//         --min-repeat-entropy 1 \
//         --min-alternate-count 2 \
//         ${intervalsOptions} \
//         ${bamTumor} \
//         ${bamNormal} > ${intervalBed.baseName}_${idSampleTumor}_vs_${idSampleNormal}.vcf
//     """
// }

// vcfFreeBayes = vcfFreeBayes.groupTuple(by:[0,1,2])

// // STEP GATK MUTECT2.1 - RAW CALLS

// process Mutect2 {
//     tag "${idSampleTumor}_vs_${idSampleNormal}-${intervalBed.baseName}"

//     label 'cpus_1'

//     input:
//         set idPatient, idSampleNormal, file(bamNormal), file(baiNormal), idSampleTumor, file(bamTumor), file(baiTumor), file(intervalBed) from pairBamMutect2
//         file(dict) from dict
//         file(fasta) from fasta
//         file(fastaFai) from fai
//         file(germlineResource) from germline_resource
//         file(germlineResourceIndex) from germline_resource_tbi
//         file(intervals) from intervals
//         file(pon) from pon
//         file(ponIndex) from pon_tbi

//     output:
//         set val("Mutect2"), idPatient, val("${idSampleTumor}_vs_${idSampleNormal}"), file("${intervalBed.baseName}_${idSampleTumor}_vs_${idSampleNormal}.vcf") into mutect2Output
//         set idPatient, idSampleNormal, idSampleTumor, file("${intervalBed.baseName}_${idSampleTumor}_vs_${idSampleNormal}.vcf.stats") optional true into intervalStatsFiles
//         set idPatient, val("${idSampleTumor}_vs_${idSampleNormal}"), file("${intervalBed.baseName}_${idSampleTumor}_vs_${idSampleNormal}.vcf.stats"), file("${intervalBed.baseName}_${idSampleTumor}_vs_${idSampleNormal}.vcf") optional true into mutect2Stats

//     when: 'mutect2' in tools

//     script:
//     // please make a panel-of-normals, using at least 40 samples
//     // https://gatkforums.broadinstitute.org/gatk/discussion/11136/how-to-call-somatic-mutations-using-gatk4-mutect2
//     PON = params.pon ? "--panel-of-normals ${pon}" : ""
//     intervalsOptions = params.no_intervals ? "" : "-L ${intervalBed}"
//     softClippedOption = params.ignore_soft_clipped_bases ? "--dont-use-soft-clipped-bases true" : ""
//     """
//     # Get raw calls
//     gatk --java-options "-Xmx${task.memory.toGiga()}g" \
//       Mutect2 \
//       -R ${fasta}\
//       -I ${bamTumor}  -tumor ${idSampleTumor} \
//       -I ${bamNormal} -normal ${idSampleNormal} \
//       ${intervalsOptions} \
//       ${softClippedOption} \
//       --germline-resource ${germlineResource} \
//       ${PON} \
//       -O ${intervalBed.baseName}_${idSampleTumor}_vs_${idSampleNormal}.vcf
//     """
// }

// mutect2Output = mutect2Output.groupTuple(by:[0,1,2])
// mutect2Stats = mutect2Stats.groupTuple(by:[0,1])

// // STEP GATK MUTECT2.2 - MERGING STATS

// process MergeMutect2Stats {
//     tag "${idSamplePair}"

//     publishDir "${params.outdir}/VariantCalling/${idSamplePair}/Mutect2", mode: params.publish_dir_mode

//     input:
//         set idPatient, idSamplePair, file(statsFiles), file(vcf) from mutect2Stats // Actual stats files and corresponding VCF chunks
//         file(dict) from dict
//         file(fasta) from fasta
//         file(fastaFai) from fai
//         file(germlineResource) from germline_resource
//         file(germlineResourceIndex) from germline_resource_tbi
//         file(intervals) from intervals

//     output:
//         set idPatient, idSamplePair, file("${idSamplePair}.vcf.gz.stats") into mergedStatsFile

//     when: 'mutect2' in tools

//     script:
//     stats = statsFiles.collect{ "-stats ${it} " }.join(' ')
//     """
//     gatk --java-options "-Xmx${task.memory.toGiga()}g" \
//         MergeMutectStats \
//         ${stats} \
//         -O ${idSamplePair}.vcf.gz.stats
//     """
// }

// // we are merging the VCFs that are called separatelly for different intervals
// // so we can have a single sorted VCF containing all the calls for a given caller

// // STEP MERGING VCF - FREEBAYES & GATK HAPLOTYPECALLER

// vcfConcatenateVCFs = vcfFreeBayes.mix(vcfFreebayesSingle, vcfGenotypeGVCFs, gvcfHaplotypeCaller)
// vcfConcatenateVCFs = vcfConcatenateVCFs.dump(tag:'VCF to merge')

// process ConcatVCF {
//     label 'cpus_8'

//     tag "${variantCaller}-${idSample}"

//     publishDir "${params.outdir}/VariantCalling/${idSample}/${"$variantCaller"}", mode: params.publish_dir_mode

//     input:
//         set variantCaller, idPatient, idSample, file(vcf) from vcfConcatenateVCFs
//         file(fastaFai) from fai
//         file(targetBED) from ch_target_bed

//     output:
//     // we have this funny *_* pattern to avoid copying the raw calls to publishdir
//         set variantCaller, idPatient, idSample, file("*_*.vcf.gz"), file("*_*.vcf.gz.tbi") into vcfConcatenated

//     when: ('haplotypecaller' in tools || 'mutect2' in tools || 'freebayes' in tools)

//     script:
//     if (variantCaller == 'HaplotypeCallerGVCF')
//         outputFile = "HaplotypeCaller_${idSample}.g.vcf"
//     else
//         outputFile = "${variantCaller}_${idSample}.vcf"
//     options = params.target_bed ? "-t ${targetBED}" : ""
//     intervalsOptions = params.no_intervals ? "-n" : ""
//     """
//     concatenateVCFs.sh -i ${fastaFai} -c ${task.cpus} -o ${outputFile} ${options} ${intervalsOptions}
//     """
// }

// vcfConcatenated = vcfConcatenated.dump(tag:'VCF')

// // STEP MERGING VCF - GATK MUTECT2 (UNFILTERED)

// mutect2Output = mutect2Output.dump(tag:'Mutect2 output VCF to merge')

// process ConcatVCF_Mutect2 {
//     label 'cpus_8'

//     tag "${idSample}"

//     publishDir "${params.outdir}/VariantCalling/${idSample}/Mutect2", mode: params.publish_dir_mode

//     input:
//         set variantCaller, idPatient, idSample, file(vcf) from mutect2Output
//         file(fastaFai) from fai
//         file(targetBED) from ch_target_bed

//     output:
//     // we have this funny *_* pattern to avoid copying the raw calls to publishdir
//         set variantCaller, idPatient, idSample, file("*_*.vcf.gz"), file("*_*.vcf.gz.tbi") into vcfConcatenatedForFilter

//     when: ('haplotypecaller' in tools || 'mutect2' in tools || 'freebayes' in tools)

//     script:
//     outputFile = "Mutect2_unfiltered_${idSample}.vcf"
//     options = params.target_bed ? "-t ${targetBED}" : ""
//     intervalsOptions = params.no_intervals ? "-n" : ""
//     """
//     concatenateVCFs.sh -i ${fastaFai} -c ${task.cpus} -o ${outputFile} ${options} ${intervalsOptions}
//     """
// }

// vcfConcatenatedForFilter = vcfConcatenatedForFilter.dump(tag:'Mutect2 unfiltered VCF')

// // STEP GATK MUTECT2.3 - GENERATING PILEUP SUMMARIES

// pairBamPileupSummaries = pairBamPileupSummaries.map{
//     idPatient, idSampleNormal, bamNormal, baiNormal, idSampleTumor, bamTumor, baiTumor, intervalBed ->
//     [idPatient, idSampleNormal, idSampleTumor, bamNormal, baiNormal, bamTumor, baiTumor, intervalBed]
// }.join(intervalStatsFiles, by:[0,1,2])

// process PileupSummariesForMutect2 {
//     tag "${idSampleTumor}_vs_${idSampleNormal}-${intervalBed.baseName}"

//     label 'cpus_1'

//     input:
//         set idPatient, idSampleNormal, idSampleTumor, file(bamNormal), file(baiNormal), file(bamTumor), file(baiTumor), file(intervalBed), file(statsFile) from pairBamPileupSummaries
//         file(germlineResource) from germline_resource
//         file(germlineResourceIndex) from germline_resource_tbi

//     output:
//         set idPatient, idSampleNormal, idSampleTumor, file("${intervalBed.baseName}_${idSampleTumor}_pileupsummaries.table") into pileupSummaries

//     when: 'mutect2' in tools

//     script:
//     intervalsOptions = params.no_intervals ? "" : "-L ${intervalBed}"
//     """
//     gatk --java-options "-Xmx${task.memory.toGiga()}g" \
//         GetPileupSummaries \
//         -I ${bamTumor} \
//         -V ${germlineResource} \
//         ${intervalsOptions} \
//         -O ${intervalBed.baseName}_${idSampleTumor}_pileupsummaries.table
//     """
// }

// pileupSummaries = pileupSummaries.groupTuple(by:[0,1,2])

// // STEP GATK MUTECT2.4 - MERGING PILEUP SUMMARIES

// process MergePileupSummaries {
//     label 'cpus_1'

//     tag "${idPatient}_${idSampleTumor}"

//     publishDir "${params.outdir}/VariantCalling/${idSampleTumor}/Mutect2", mode: params.publish_dir_mode

//     input:
//         set idPatient, idSampleNormal, idSampleTumor, file(pileupSums) from pileupSummaries
//         file(dict) from dict

//     output:
//         set idPatient, idSampleNormal, idSampleTumor, file("${idSampleTumor}_pileupsummaries.table") into mergedPileupFile

//     when: 'mutect2' in tools

//     script:
//     allPileups = pileupSums.collect{ "-I ${it} " }.join(' ')
//     """
//     gatk --java-options "-Xmx${task.memory.toGiga()}g" \
//         GatherPileupSummaries \
//         --sequence-dictionary ${dict} \
//         ${allPileups} \
//         -O ${idSampleTumor}_pileupsummaries.table
//     """
// }

// // STEP GATK MUTECT2.5 - CALCULATING CONTAMINATION

// pairBamCalculateContamination = pairBamCalculateContamination.map{
//     idPatient, idSampleNormal, bamNormal, baiNormal, idSampleTumor, bamTumor, baiTumor ->
//     [idPatient, idSampleNormal, idSampleTumor, bamNormal, baiNormal, bamTumor, baiTumor]
// }.join(mergedPileupFile, by:[0,1,2])

// process CalculateContamination {
//     label 'cpus_1'

//     tag "${idSampleTumor}_vs_${idSampleNormal}"

//     publishDir "${params.outdir}/VariantCalling/${idSampleTumor}/Mutect2", mode: params.publish_dir_mode

//     input:
//         set idPatient, idSampleNormal, idSampleTumor, file(bamNormal), file(baiNormal), file(bamTumor), file(baiTumor), file(mergedPileup) from pairBamCalculateContamination

//      output:
//         set idPatient, val("${idSampleTumor}_vs_${idSampleNormal}"), file("${idSampleTumor}_contamination.table") into contaminationTable

//     when: 'mutect2' in tools

//     script:
//     """
//     # calculate contamination
//     gatk --java-options "-Xmx${task.memory.toGiga()}g" \
//         CalculateContamination \
//         -I ${idSampleTumor}_pileupsummaries.table \
//         -O ${idSampleTumor}_contamination.table
//     """
// }

// // STEP GATK MUTECT2.6 - FILTERING CALLS

// mutect2CallsToFilter = vcfConcatenatedForFilter.map{
//     variantCaller, idPatient, idSamplePair, vcf, tbi ->
//     [idPatient, idSamplePair, vcf, tbi]
// }.join(mergedStatsFile, by:[0,1]).join(contaminationTable, by:[0,1])

// process FilterMutect2Calls {
//     label 'cpus_1'

//     tag "${idSamplePair}"

//     publishDir "${params.outdir}/VariantCalling/${idSamplePair}/Mutect2", mode: params.publish_dir_mode

//     input:
//         set idPatient, idSamplePair, file(unfiltered), file(unfilteredIndex), file(stats), file(contaminationTable) from mutect2CallsToFilter
//         file(dict) from dict
//         file(fasta) from fasta
//         file(fastaFai) from fai
//         file(germlineResource) from germline_resource
//         file(germlineResourceIndex) from germline_resource_tbi
//         file(intervals) from intervals

//     output:
//         set val("Mutect2"), idPatient, idSamplePair, file("Mutect2_filtered_${idSamplePair}.vcf.gz"), file("Mutect2_filtered_${idSamplePair}.vcf.gz.tbi"), file("Mutect2_filtered_${idSamplePair}.vcf.gz.filteringStats.tsv") into filteredMutect2Output

//     when: 'mutect2' in tools

//     script:
//     """
//     # do the actual filtering
//     gatk --java-options "-Xmx${task.memory.toGiga()}g" \
//         FilterMutectCalls \
//         -V ${unfiltered} \
//         --contamination-table ${contaminationTable} \
//         --stats ${stats} \
//         -R ${fasta} \
//         -O Mutect2_filtered_${idSamplePair}.vcf.gz
//     """
// }

// // STEP SENTIEON TNSCOPE

// process Sentieon_TNscope {
//     label 'cpus_max'
//     label 'memory_max'
//     label 'sentieon'

//     tag "${idSampleTumor}_vs_${idSampleNormal}"

//     input:
//         set idPatient, idSampleNormal, file(bamNormal), file(baiNormal), file(recalNormal), idSampleTumor, file(bamTumor), file(baiTumor), file(recalTumor) from bam_pair_sentieon_TNscope
//         file(dict) from dict
//         file(fasta) from fasta
//         file(fastaFai) from fai
//         file(dbsnp) from dbsnp
//         file(dbsnpIndex) from dbsnp_tbi
//         file(pon) from pon
//         file(ponIndex) from pon_tbi

//     output:
//         set val("SentieonTNscope"), idPatient, val("${idSampleTumor}_vs_${idSampleNormal}"), file("*.vcf") into vcf_sentieon_TNscope

//     when: 'tnscope' in tools && params.sentieon

//     script:
//     PON = params.pon ? "--pon ${pon}" : ""
//     """
//     sentieon driver \
//         -t ${task.cpus} \
//         -r ${fasta} \
//         -i ${bamTumor} \
//         -q ${recalTumor} \
//         -i ${bamNormal} \
//         -q ${recalNormal} \
//         --algo TNscope \
//         --tumor_sample ${idSampleTumor} \
//         --normal_sample ${idSampleNormal} \
//         --dbsnp ${dbsnp} \
//         ${PON} \
//         TNscope_${idSampleTumor}_vs_${idSampleNormal}.vcf
//     """
// }

// vcf_sentieon_TNscope = vcf_sentieon_TNscope.dump(tag:'Sentieon TNscope')

// vcf_sentieon = vcf_sentieon_DNAseq.mix(vcf_sentieon_DNAscope, vcf_sentieon_DNAscope_SV, vcf_sentieon_TNscope)

// process CompressSentieonVCF {
//     tag "${idSample} - ${vcf}"

//     publishDir "${params.outdir}/VariantCalling/${idSample}/${variantCaller}", mode: params.publish_dir_mode

//     input:
//         set variantCaller, idPatient, idSample, file(vcf) from vcf_sentieon

//     output:
//         set variantCaller, idPatient, idSample, file("*.vcf.gz"), file("*.vcf.gz.tbi") into vcf_sentieon_compressed

//     script:
//     """
//     bgzip < ${vcf} > ${vcf}.gz
//     tabix ${vcf}.gz
//     """
// }

// vcf_sentieon_compressed = vcf_sentieon_compressed.dump(tag:'Sentieon VCF indexed')

// // STEP STRELKA.2 - SOMATIC PAIR

// process Strelka {
//     label 'cpus_max'
//     label 'memory_max'

//     tag "${idSampleTumor}_vs_${idSampleNormal}"

//     publishDir "${params.outdir}/VariantCalling/${idSampleTumor}_vs_${idSampleNormal}/Strelka", mode: params.publish_dir_mode

//     input:
//         set idPatient, idSampleNormal, file(bamNormal), file(baiNormal), idSampleTumor, file(bamTumor), file(baiTumor) from pairBamStrelka
//         file(dict) from dict
//         file(fasta) from fasta
//         file(fastaFai) from fai
//         file(targetBED) from ch_target_bed

//     output:
//         set val("Strelka"), idPatient, val("${idSampleTumor}_vs_${idSampleNormal}"), file("*.vcf.gz"), file("*.vcf.gz.tbi") into vcfStrelka

//     when: 'strelka' in tools

//     script:
//     beforeScript = params.target_bed ? "bgzip --threads ${task.cpus} -c ${targetBED} > call_targets.bed.gz ; tabix call_targets.bed.gz" : ""
//     options = params.target_bed ? "--exome --callRegions call_targets.bed.gz" : ""
//     """
//     ${beforeScript}
//     configureStrelkaSomaticWorkflow.py \
//         --tumor ${bamTumor} \
//         --normal ${bamNormal} \
//         --referenceFasta ${fasta} \
//         ${options} \
//         --runDir Strelka

//     python Strelka/runWorkflow.py -m local -j ${task.cpus}

//     mv Strelka/results/variants/somatic.indels.vcf.gz \
//         Strelka_${idSampleTumor}_vs_${idSampleNormal}_somatic_indels.vcf.gz
//     mv Strelka/results/variants/somatic.indels.vcf.gz.tbi \
//         Strelka_${idSampleTumor}_vs_${idSampleNormal}_somatic_indels.vcf.gz.tbi
//     mv Strelka/results/variants/somatic.snvs.vcf.gz \
//         Strelka_${idSampleTumor}_vs_${idSampleNormal}_somatic_snvs.vcf.gz
//     mv Strelka/results/variants/somatic.snvs.vcf.gz.tbi \
//         Strelka_${idSampleTumor}_vs_${idSampleNormal}_somatic_snvs.vcf.gz.tbi
//     """
// }

// vcfStrelka = vcfStrelka.dump(tag:'Strelka')

// // STEP MANTA.2 - SOMATIC PAIR

// process Manta {
//     label 'cpus_max'
//     label 'memory_max'

//     tag "${idSampleTumor}_vs_${idSampleNormal}"

//     publishDir "${params.outdir}/VariantCalling/${idSampleTumor}_vs_${idSampleNormal}/Manta", mode: params.publish_dir_mode

//     input:
//         set idPatient, idSampleNormal, file(bamNormal), file(baiNormal), idSampleTumor, file(bamTumor), file(baiTumor) from pairBamManta
//         file(fasta) from fasta
//         file(fastaFai) from fai
//         file(targetBED) from ch_target_bed

//     output:
//         set val("Manta"), idPatient, val("${idSampleTumor}_vs_${idSampleNormal}"), file("*.vcf.gz"), file("*.vcf.gz.tbi") into vcfManta
//         set idPatient, idSampleNormal, idSampleTumor, file("*.candidateSmallIndels.vcf.gz"), file("*.candidateSmallIndels.vcf.gz.tbi") into mantaToStrelka

//     when: 'manta' in tools

//     script:
//     beforeScript = params.target_bed ? "bgzip --threads ${task.cpus} -c ${targetBED} > call_targets.bed.gz ; tabix call_targets.bed.gz" : ""
//     options = params.target_bed ? "--exome --callRegions call_targets.bed.gz" : ""
//     """
//     ${beforeScript}
//     configManta.py \
//         --normalBam ${bamNormal} \
//         --tumorBam ${bamTumor} \
//         --reference ${fasta} \
//         ${options} \
//         --runDir Manta

//     python Manta/runWorkflow.py -m local -j ${task.cpus}

//     mv Manta/results/variants/candidateSmallIndels.vcf.gz \
//         Manta_${idSampleTumor}_vs_${idSampleNormal}.candidateSmallIndels.vcf.gz
//     mv Manta/results/variants/candidateSmallIndels.vcf.gz.tbi \
//         Manta_${idSampleTumor}_vs_${idSampleNormal}.candidateSmallIndels.vcf.gz.tbi
//     mv Manta/results/variants/candidateSV.vcf.gz \
//         Manta_${idSampleTumor}_vs_${idSampleNormal}.candidateSV.vcf.gz
//     mv Manta/results/variants/candidateSV.vcf.gz.tbi \
//         Manta_${idSampleTumor}_vs_${idSampleNormal}.candidateSV.vcf.gz.tbi
//     mv Manta/results/variants/diploidSV.vcf.gz \
//         Manta_${idSampleTumor}_vs_${idSampleNormal}.diploidSV.vcf.gz
//     mv Manta/results/variants/diploidSV.vcf.gz.tbi \
//         Manta_${idSampleTumor}_vs_${idSampleNormal}.diploidSV.vcf.gz.tbi
//     mv Manta/results/variants/somaticSV.vcf.gz \
//         Manta_${idSampleTumor}_vs_${idSampleNormal}.somaticSV.vcf.gz
//     mv Manta/results/variants/somaticSV.vcf.gz.tbi \
//         Manta_${idSampleTumor}_vs_${idSampleNormal}.somaticSV.vcf.gz.tbi
//     """
// }

// vcfManta = vcfManta.dump(tag:'Manta')

// // Remmaping channels to match input for StrelkaBP
// pairBamStrelkaBP = pairBamStrelkaBP.map {
//     idPatientNormal, idSampleNormal, bamNormal, baiNormal, idSampleTumor, bamTumor, baiTumor ->
//     [idPatientNormal, idSampleNormal, idSampleTumor, bamNormal, baiNormal, bamTumor, baiTumor]
// }.join(mantaToStrelka, by:[0,1,2]).map {
//     idPatientNormal, idSampleNormal, idSampleTumor, bamNormal, baiNormal, bamTumor, baiTumor, mantaCSI, mantaCSIi ->
//     [idPatientNormal, idSampleNormal, bamNormal, baiNormal, idSampleTumor, bamTumor, baiTumor, mantaCSI, mantaCSIi]
// }

// // STEP STRELKA.3 - SOMATIC PAIR - BEST PRACTICES

// process StrelkaBP {
//     label 'cpus_max'
//     label 'memory_max'

//     tag "${idSampleTumor}_vs_${idSampleNormal}"

//     publishDir "${params.outdir}/VariantCalling/${idSampleTumor}_vs_${idSampleNormal}/Strelka", mode: params.publish_dir_mode

//     input:
//         set idPatient, idSampleNormal, file(bamNormal), file(baiNormal), idSampleTumor, file(bamTumor), file(baiTumor), file(mantaCSI), file(mantaCSIi) from pairBamStrelkaBP
//         file(dict) from dict
//         file(fasta) from fasta
//         file(fastaFai) from fai
//         file(targetBED) from ch_target_bed

//     output:
//         set val("Strelka"), idPatient, val("${idSampleTumor}_vs_${idSampleNormal}"), file("*.vcf.gz"), file("*.vcf.gz.tbi") into vcfStrelkaBP

//     when: 'strelka' in tools && 'manta' in tools && !params.no_strelka_bp

//     script:
//     beforeScript = params.target_bed ? "bgzip --threads ${task.cpus} -c ${targetBED} > call_targets.bed.gz ; tabix call_targets.bed.gz" : ""
//     options = params.target_bed ? "--exome --callRegions call_targets.bed.gz" : ""
//     """
//     ${beforeScript}
//     configureStrelkaSomaticWorkflow.py \
//         --tumor ${bamTumor} \
//         --normal ${bamNormal} \
//         --referenceFasta ${fasta} \
//         --indelCandidates ${mantaCSI} \
//         ${options} \
//         --runDir Strelka

//     python Strelka/runWorkflow.py -m local -j ${task.cpus}

//     mv Strelka/results/variants/somatic.indels.vcf.gz \
//         StrelkaBP_${idSampleTumor}_vs_${idSampleNormal}_somatic_indels.vcf.gz
//     mv Strelka/results/variants/somatic.indels.vcf.gz.tbi \
//         StrelkaBP_${idSampleTumor}_vs_${idSampleNormal}_somatic_indels.vcf.gz.tbi
//     mv Strelka/results/variants/somatic.snvs.vcf.gz \
//         StrelkaBP_${idSampleTumor}_vs_${idSampleNormal}_somatic_snvs.vcf.gz
//     mv Strelka/results/variants/somatic.snvs.vcf.gz.tbi \
//         StrelkaBP_${idSampleTumor}_vs_${idSampleNormal}_somatic_snvs.vcf.gz.tbi
//     """
// }

// vcfStrelkaBP = vcfStrelkaBP.dump(tag:'Strelka BP')

// // STEP CNVkit

// process CNVkit {
//     tag "${idSampleTumor}_vs_${idSampleNormal}"

//     publishDir "${params.outdir}/VariantCalling/${idSampleTumor}_vs_${idSampleNormal}/CNVkit", mode: params.publish_dir_mode

//     input:
//         set idPatient, idSampleNormal, file(bamNormal), file(baiNormal), idSampleTumor, file(bamTumor), file(baiTumor) from pairBamCNVkit
//         file(targetBED) from ch_target_bed
//         file(fasta) from fasta

//     output:
//         set idPatient, idSampleNormal, idSampleTumor, file("${idSampleTumor}*"), file("${idSampleNormal}*") into cnvkitOut

//     when: 'cnvkit' in tools && params.target_bed

//     script:
//     """
//     cnvkit.py \
//       batch \
//       ${bamTumor} \
//       --normal ${bamNormal} \
//       --targets ${targetBED} \
//       --fasta ${fasta} \
//       --output-reference output_reference.cnn \
//       --output-dir ./ \
//       --diagram \
//       --scatter
//     """
// }

// // STEP MSISENSOR.1 - SCAN

// // Scan reference genome for microsatellites
// process MSIsensor_scan {
//     label 'cpus_1'
//     label 'memory_max'

//     tag "${fasta}"

//     input:
//     file(fasta) from fasta
//     file(fastaFai) from fai

//     output:
//     file "microsatellites.list" into msi_scan_ch

//     when: 'msisensor' in tools

//     script:
//     """
//     msisensor scan -d ${fasta} -o microsatellites.list
//     """
// }

// // STEP MSISENSOR.2 - SCORE

// // Score the normal vs somatic pair of bams

// process MSIsensor_msi {
//     label 'cpus_4'
//     label 'memory_max'

//     tag "${idSampleTumor}_vs_${idSampleNormal}"

//     publishDir "${params.outdir}/VariantCalling/${idSampleTumor}_vs_${idSampleNormal}/MSIsensor", mode: params.publish_dir_mode

//     input:
//         set idPatient, idSampleNormal, file(bamNormal), file(baiNormal), idSampleTumor, file(bamTumor), file(baiTumor) from pairBamMsisensor
//         file msiSites from msi_scan_ch

//     output:
//         set val("Msisensor"), idPatient, file("${idSampleTumor}_vs_${idSampleNormal}_msisensor"), file("${idSampleTumor}_vs_${idSampleNormal}_msisensor_dis"), file("${idSampleTumor}_vs_${idSampleNormal}_msisensor_germline"), file("${idSampleTumor}_vs_${idSampleNormal}_msisensor_somatic") into msisensor_out_ch

//     when: 'msisensor' in tools

//     script:
//     """
//     msisensor msi -d ${msiSites} \
//                   -b 4 \
//                   -n ${bamNormal} \
//                   -t ${bamTumor} \
//                   -o ${idSampleTumor}_vs_${idSampleNormal}_msisensor
//     """
// }

// // STEP ASCAT.1 - ALLELECOUNTER

// // Run commands and code from Malin Larsson
// // Based on Jesper Eisfeldt's code
// process AlleleCounter {
//     label 'memory_singleCPU_2_task'

//     tag "${idSample}"

//     input:
//         set idPatient, idSample, file(bam), file(bai) from bamAscat
//         file(acLoci) from ch_ac_loci
//         file(dict) from dict
//         file(fasta) from fasta
//         file(fastaFai) from fai

//     output:
//         set idPatient, idSample, file("${idSample}.alleleCount") into alleleCounterOut

//     when: 'ascat' in tools

//     script:
//     """
//     alleleCounter \
//         -l ${acLoci} \
//         -r ${fasta} \
//         -b ${bam} \
//         -o ${idSample}.alleleCount;
//     """
// }

// alleleCountOutNormal = Channel.create()
// alleleCountOutTumor = Channel.create()

// alleleCounterOut
//     .choice(alleleCountOutTumor, alleleCountOutNormal) {status_map[it[0], it[1]] == 0 ? 1 : 0}

// alleleCounterOut = alleleCountOutNormal.combine(alleleCountOutTumor, by:0)

// alleleCounterOut = alleleCounterOut.map {
//     idPatientNormal, idSampleNormal, alleleCountOutNormal,
//     idSampleTumor, alleleCountOutTumor ->
//     [idPatientNormal, idSampleNormal, idSampleTumor, alleleCountOutNormal, alleleCountOutTumor]
// }

// // STEP ASCAT.2 - CONVERTALLELECOUNTS

// // R script from Malin Larssons bitbucket repo:
// // https://bitbucket.org/malinlarsson/somatic_wgs_pipeline
// process ConvertAlleleCounts {
//     label 'memory_singleCPU_2_task'

//     tag "${idSampleTumor}_vs_${idSampleNormal}"

//     publishDir "${params.outdir}/VariantCalling/${idSampleTumor}_vs_${idSampleNormal}/ASCAT", mode: params.publish_dir_mode

//     input:
//         set idPatient, idSampleNormal, idSampleTumor, file(alleleCountNormal), file(alleleCountTumor) from alleleCounterOut

//     output:
//         set idPatient, idSampleNormal, idSampleTumor, file("${idSampleNormal}.BAF"), file("${idSampleNormal}.LogR"), file("${idSampleTumor}.BAF"), file("${idSampleTumor}.LogR") into convertAlleleCountsOut

//     when: 'ascat' in tools

//     script:
//     gender = gender_map[idPatient]
//     """
//     convertAlleleCounts.r ${idSampleTumor} ${alleleCountTumor} ${idSampleNormal} ${alleleCountNormal} ${gender}
//     """
// }

// // STEP ASCAT.3 - ASCAT

// // R scripts from Malin Larssons bitbucket repo:
// // https://bitbucket.org/malinlarsson/somatic_wgs_pipeline
// process Ascat {
//     label 'memory_singleCPU_2_task'

//     tag "${idSampleTumor}_vs_${idSampleNormal}"

//     publishDir "${params.outdir}/VariantCalling/${idSampleTumor}_vs_${idSampleNormal}/ASCAT", mode: params.publish_dir_mode

//     input:
//         set idPatient, idSampleNormal, idSampleTumor, file(bafNormal), file(logrNormal), file(bafTumor), file(logrTumor) from convertAlleleCountsOut
//         file(acLociGC) from ch_ac_loci_gc

//     output:
//         set val("ASCAT"), idPatient, idSampleNormal, idSampleTumor, file("${idSampleTumor}.*.{png,txt}") into ascatOut

//     when: 'ascat' in tools

//     script:
//     gender = gender_map[idPatient]
//     purity_ploidy = (params.ascat_purity && params.ascat_ploidy) ? "--purity ${params.ascat_purity} --ploidy ${params.ascat_ploidy}" : ""
//     """
//     for f in *BAF *LogR; do sed 's/chr//g' \$f > tmpFile; mv tmpFile \$f;done
//     run_ascat.r \
//         --tumorbaf ${bafTumor} \
//         --tumorlogr ${logrTumor} \
//         --normalbaf ${bafNormal} \
//         --normallogr ${logrNormal} \
//         --tumorname ${idSampleTumor} \
//         --basedir ${baseDir} \
//         --gcfile ${acLociGC} \
//         --gender ${gender} \
//         ${purity_ploidy}
//     """
// }

// ascatOut.dump(tag:'ASCAT')

// // STEP MPILEUP.1

// process Mpileup {
//     label 'cpus_1'
//     label 'memory_singleCPU_2_task'

//     tag "${idSample}-${intervalBed.baseName}"

//     publishDir params.outdir, mode: params.publish_dir_mode, saveAs: { it == "${idSample}.pileup" ? "VariantCalling/${idSample}/Control-FREEC/${it}" : null }

//     input:
//         set idPatient, idSample, file(bam), file(bai), file(intervalBed) from bamMpileup
//         file(fasta) from fasta
//         file(fastaFai) from fai

//     output:
//         set idPatient, idSample, file("${prefix}${idSample}.pileup") into mpileupMerge
//         set idPatient, idSample into tsv_mpileup

//     when: 'controlfreec' in tools || 'mpileup' in tools

//     script:
//     prefix = params.no_intervals ? "" : "${intervalBed.baseName}_"
//     intervalsOptions = params.no_intervals ? "" : "-l ${intervalBed}"

//     """
//     # Control-FREEC reads uncompresses the zipped file TWICE in single-threaded mode.
//     # we are therefore not using compressed pileups here
//     samtools mpileup \
//         -f ${fasta} ${bam} \
//         ${intervalsOptions} > ${prefix}${idSample}.pileup
//     """
// }

// (tsv_mpileup, tsv_mpileup_sample) = tsv_mpileup.groupTuple(by:[0, 1]).into(2)

// // Creating a TSV file to restart from this step
// tsv_mpileup.map { idPatient, idSample ->
//     gender = gender_map[idPatient]
//     status = status_map[idPatient, idSample]
//     mpileup = "${params.outdir}/VariantCalling/${idSample}/Control-FREEC/${idSample}.pileup"
//     "${idPatient}\t${gender}\t${status}\t${idSample}\t${mpileup}\n"
// }.collectFile(
//     name: 'control-freec_mpileup.tsv', sort: true, storeDir: "${params.outdir}/VariantCalling/TSV"
// )

// tsv_mpileup_sample
//     .collectFile(storeDir: "${params.outdir}/VariantCalling/TSV") {
//         idPatient, idSample ->
//         status = status_map[idPatient, idSample]
//         gender = gender_map[idPatient]
//         mpileup = "${params.outdir}/VariantCalling/${idSample}/Control-FREEC/${idSample}.pileup"
//         ["control-freec_mpileup_${idSample}.tsv", "${idPatient}\t${gender}\t${status}\t${idSample}\t${mpileup}\n"]
// }

// if (!params.no_intervals) {
//     mpileupMerge = mpileupMerge.groupTuple(by:[0, 1])
//     mpileupNoInt = Channel.empty()
// } else {
//     (mpileupMerge, mpileupNoInt) = mpileupMerge.into(2)
//     mpileupMerge.close()
// }

// // STEP MPILEUP.2 - MERGE
// process MergeMpileup {
//     label 'cpus_1'

//     tag "${idSample}"

//     publishDir params.outdir, mode: params.publish_dir_mode, saveAs: { it == "${idSample}.pileup" ? "VariantCalling/${idSample}/Control-FREEC/${it}" : null }

//     input:
//         set idPatient, idSample, file(mpileup) from mpileupMerge

//     output:
//         set idPatient, idSample, file("${idSample}.pileup") into mpileupOut

//     when: !(params.no_intervals) && 'controlfreec' in tools || 'mpileup' in tools

//     script:
//     """
//     for i in `ls -1v *.pileup`;
//         do cat \$i >> ${idSample}.pileup
//     done
//     """
// }

// mpileupOut = mpileupOut.mix(mpileupNoInt)
// mpileupOut = mpileupOut.dump(tag:'mpileup')

// mpileupOutNormal = Channel.create()
// mpileupOutTumor = Channel.create()

// if (step == 'controlfreec') mpileupOut = input_sample

// mpileupOut
//     .choice(mpileupOutTumor, mpileupOutNormal) {status_map[it[0], it[1]] == 0 ? 1 : 0}

// mpileupOut = mpileupOutNormal.combine(mpileupOutTumor, by:0)

// mpileupOut = mpileupOut.map {
//     idPatientNormal, idSampleNormal, mpileupOutNormal,
//     idSampleTumor, mpileupOutTumor ->
//     [idPatientNormal, idSampleNormal, idSampleTumor, mpileupOutNormal, mpileupOutTumor]
// }

// // STEP CONTROLFREEC.1 - CONTROLFREEC

// process ControlFREEC {
//     label 'cpus_max'
//     //label 'memory_singleCPU_2_task'

//     tag "${idSampleTumor}_vs_${idSampleNormal}"

//     publishDir "${params.outdir}/VariantCalling/${idSampleTumor}_vs_${idSampleNormal}/Control-FREEC", mode: params.publish_dir_mode

//     input:
//         set idPatient, idSampleNormal, idSampleTumor, file(mpileupNormal), file(mpileupTumor) from mpileupOut
//         file(chrDir) from ch_chr_dir
//         file(mappability) from ch_mappability
//         file(chrLength) from ch_chr_length
//         file(dbsnp) from dbsnp
//         file(dbsnpIndex) from dbsnp_tbi
//         file(fasta) from fasta
//         file(fastaFai) from fai

//     output:
//         set idPatient, idSampleNormal, idSampleTumor, file("${idSampleTumor}.pileup_CNVs"), file("${idSampleTumor}.pileup_ratio.txt"), file("${idSampleTumor}.pileup_normal_CNVs"), file("${idSampleTumor}.pileup_normal_ratio.txt"), file("${idSampleTumor}.pileup_BAF.txt"), file("${idSampleNormal}.pileup_BAF.txt") into controlFreecViz
//         set file("*.pileup*"), file("${idSampleTumor}_vs_${idSampleNormal}.config.txt") into controlFreecOut

//     when: 'controlfreec' in tools

//     script:
//     config = "${idSampleTumor}_vs_${idSampleNormal}.config.txt"
//     gender = gender_map[idPatient]
//     // if we are using coefficientOfVariation, we must delete the window parameter
//     // it is "window = 20000" in the default settings, without coefficientOfVariation set,
//     // but we do not like it. Note, it is not written in stone
//     coeff_or_window = params.cf_window ? "window = ${params.cf_window}" : "coefficientOfVariation = ${params.cf_coeff}"

//     """
//     touch ${config}
//     echo "[general]" >> ${config}
//     echo "BedGraphOutput = TRUE" >> ${config}
//     echo "chrFiles = \${PWD}/${chrDir.fileName}" >> ${config}
//     echo "chrLenFile = \${PWD}/${chrLength.fileName}" >> ${config}
//     echo "gemMappabilityFile = \${PWD}/${mappability}" >> ${config}
//     echo "${coeff_or_window}" >> ${config}
//     echo "contaminationAdjustment = TRUE" >> ${config}
//     echo "forceGCcontentNormalization = 1" >> ${config}
//     echo "maxThreads = ${task.cpus}" >> ${config}
//     echo "minimalSubclonePresence = 20" >> ${config}
//     echo "ploidy = ${params.cf_ploidy}" >> ${config}
//     echo "sex = ${gender}" >> ${config}
//     echo "" >> ${config}

//     echo "[control]" >> ${config}
//     echo "inputFormat = pileup" >> ${config}
//     echo "mateFile = \${PWD}/${mpileupNormal}" >> ${config}
//     echo "mateOrientation = FR" >> ${config}
//     echo "" >> ${config}

//     echo "[sample]" >> ${config}
//     echo "inputFormat = pileup" >> ${config}
//     echo "mateFile = \${PWD}/${mpileupTumor}" >> ${config}
//     echo "mateOrientation = FR" >> ${config}
//     echo "" >> ${config}

//     echo "[BAF]" >> ${config}
//     echo "SNPfile = ${dbsnp.fileName}" >> ${config}

//     freec -conf ${config}
//     """
// }

// controlFreecOut.dump(tag:'ControlFREEC')

// // STEP CONTROLFREEC.3 - VISUALIZATION

// process ControlFreecViz {
//     label 'memory_singleCPU_2_task'

//     tag "${idSampleTumor}_vs_${idSampleNormal}"

//     publishDir "${params.outdir}/VariantCalling/${idSampleTumor}_vs_${idSampleNormal}/Control-FREEC", mode: params.publish_dir_mode

//     input:
//         set idPatient, idSampleNormal, idSampleTumor, file(cnvTumor), file(ratioTumor), file(cnvNormal), file(ratioNormal), file(bafTumor), file(bafNormal) from controlFreecViz

//     output:
//         set file("*.txt"), file("*.png"), file("*.bed") into controlFreecVizOut

//     when: 'controlfreec' in tools

//     """
//     echo "Shaping CNV files to make sure we can assess significance"
//     awk 'NF==9{print}' ${cnvTumor} > TUMOR.CNVs
//     awk 'NF==7{print}' ${cnvNormal} > NORMAL.CNVs

//     echo "############### Calculating significance values for TUMOR CNVs #############"
//     cat /opt/conda/envs/nf-core-sarek-${workflow.manifest.version}/bin/assess_significance.R | R --slave --args TUMOR.CNVs ${ratioTumor}

//     echo "############### Calculating significance values for NORMAL CNVs ############"
//     cat /opt/conda/envs/nf-core-sarek-${workflow.manifest.version}/bin/assess_significance.R | R --slave --args NORMAL.CNVs ${ratioNormal}

//     echo "############### Creating graph for TUMOR ratios ###############"
//     cat /opt/conda/envs/nf-core-sarek-${workflow.manifest.version}/bin/makeGraph.R | R --slave --args 2 ${ratioTumor} ${bafTumor}

//     echo "############### Creating graph for NORMAL ratios ##############"
//     cat /opt/conda/envs/nf-core-sarek-${workflow.manifest.version}/bin/makeGraph.R | R --slave --args 2 ${ratioNormal} ${bafNormal}

//     echo "############### Creating BED files for TUMOR ##############"
//     perl /opt/conda/envs/nf-core-sarek-${workflow.manifest.version}/bin/freec2bed.pl -f ${ratioTumor} > ${idSampleTumor}.bed

//     echo "############### Creating BED files for NORMAL #############"
//     perl /opt/conda/envs/nf-core-sarek-${workflow.manifest.version}/bin/freec2bed.pl -f ${ratioNormal} > ${idSampleNormal}.bed
//     """
// }

// controlFreecVizOut.dump(tag:'ControlFreecViz')

// // Remapping channels for QC and annotation

// (vcfStrelkaIndels, vcfStrelkaSNVS) = vcfStrelka.into(2)
// (vcfStrelkaBPIndels, vcfStrelkaBPSNVS) = vcfStrelkaBP.into(2)
// (vcfMantaSomaticSV, vcfMantaDiploidSV) = vcfManta.into(2)

// vcfKeep = Channel.empty().mix(
//     filteredMutect2Output.map{
//         variantCaller, idPatient, idSample, vcf, tbi, tsv ->
//         [variantcaller, idSample, vcf]
//     },
//     vcfConcatenated.map{
//         variantcaller, idPatient, idSample, vcf, tbi ->
//         [variantcaller, idSample, vcf]
//     },
//     vcf_sentieon_compressed.map {
//         variantcaller, idPatient, idSample, vcf, tbi ->
//         [variantcaller, idSample, vcf]
//     },
//     vcfStrelkaSingle.map {
//         variantcaller, idPatient, idSample, vcf, tbi ->
//         [variantcaller, idSample, vcf[1]]
//     },
//     vcfMantaSingle.map {
//         variantcaller, idPatient, idSample, vcf, tbi ->
//         [variantcaller, idSample, vcf[2]]
//     },
//     vcfMantaDiploidSV.map {
//         variantcaller, idPatient, idSample, vcf, tbi ->
//         [variantcaller, idSample, vcf[2]]
//     },
//     vcfMantaSomaticSV.map {
//         variantcaller, idPatient, idSample, vcf, tbi ->
//         [variantcaller, idSample, vcf[3]]
//     },
//     vcfStrelkaIndels.map {
//         variantcaller, idPatient, idSample, vcf, tbi ->
//         [variantcaller, idSample, vcf[0]]
//     },
//     vcfStrelkaSNVS.map {
//         variantcaller, idPatient, idSample, vcf, tbi ->
//         [variantcaller, idSample, vcf[1]]
//     },
//     vcfStrelkaBPIndels.map {
//         variantcaller, idPatient, idSample, vcf, tbi ->
//         [variantcaller, idSample, vcf[0]]
//     },
//     vcfStrelkaBPSNVS.map {
//         variantcaller, idPatient, idSample, vcf, tbi ->
//         [variantcaller, idSample, vcf[1]]
//     },
//     vcfTIDDIT.map {
//         variantcaller, idPatient, idSample, vcf, tbi ->
//         [variantcaller, idSample, vcf]
//     })

// (vcfBCFtools, vcfVCFtools, vcfAnnotation) = vcfKeep.into(3)

// // STEP VCF.QC

// process BcftoolsStats {
//     label 'cpus_1'

//     tag "${variantCaller} - ${vcf}"

//     publishDir "${params.outdir}/Reports/${idSample}/BCFToolsStats", mode: params.publish_dir_mode

//     input:
//         set variantCaller, idSample, file(vcf) from vcfBCFtools

//     output:
//         file ("*.bcf.tools.stats.out") into bcftoolsReport

//     when: !('bcftools' in skipQC)

//     script:
//     """
//     bcftools stats ${vcf} > ${reduceVCF(vcf.fileName)}.bcf.tools.stats.out
//     """
// }

// bcftoolsReport = bcftoolsReport.dump(tag:'BCFTools')

// process Vcftools {
//     label 'cpus_1'

//     tag "${variantCaller} - ${vcf}"

//     publishDir "${params.outdir}/Reports/${idSample}/VCFTools", mode: params.publish_dir_mode

//     input:
//         set variantCaller, idSample, file(vcf) from vcfVCFtools

//     output:
//         file ("${reduceVCF(vcf.fileName)}.*") into vcftoolsReport

//     when: !('vcftools' in skipQC)

//     script:
//     """
//     vcftools \
//     --gzvcf ${vcf} \
//     --TsTv-by-count \
//     --out ${reduceVCF(vcf.fileName)}

//     vcftools \
//     --gzvcf ${vcf} \
//     --TsTv-by-qual \
//     --out ${reduceVCF(vcf.fileName)}

//     vcftools \
//     --gzvcf ${vcf} \
//     --FILTER-summary \
//     --out ${reduceVCF(vcf.fileName)}
//     """
// }

// vcftoolsReport = vcftoolsReport.dump(tag:'VCFTools')

// /*
// ================================================================================
//                                    ANNOTATION
// ================================================================================
// */

// if (step == 'annotate') {
//     vcfToAnnotate = Channel.create()
//     vcfNoAnnotate = Channel.create()

//     if (tsv_path == []) {
//     // Sarek, by default, annotates all available vcfs that it can find in the VariantCalling directory
//     // Excluding vcfs from FreeBayes, and g.vcf from HaplotypeCaller
//     // Basically it's: results/VariantCalling/*/{HaplotypeCaller,Manta,Mutect2,SentieonDNAseq,SentieonDNAscope,SentieonTNscope,Strelka,TIDDIT}/*.vcf.gz
//     // Without *SmallIndels.vcf.gz from Manta, and *.genome.vcf.gz from Strelka
//     // The small snippet `vcf.minus(vcf.fileName)[-2]` catches idSample
//     // This field is used to output final annotated VCFs in the correct directory
//       Channel.empty().mix(
//         Channel.fromPath("${params.outdir}/VariantCalling/*/HaplotypeCaller/*.vcf.gz")
//           .flatten().map{vcf -> ['HaplotypeCaller', vcf.minus(vcf.fileName)[-2].toString(), vcf]},
//         Channel.fromPath("${params.outdir}/VariantCalling/*/Manta/*[!candidate]SV.vcf.gz")
//           .flatten().map{vcf -> ['Manta', vcf.minus(vcf.fileName)[-2].toString(), vcf]},
//         Channel.fromPath("${params.outdir}/VariantCalling/*/Mutect2/*.vcf.gz")
//           .flatten().map{vcf -> ['Mutect2', vcf.minus(vcf.fileName)[-2].toString(), vcf]},
//         Channel.fromPath("${params.outdir}/VariantCalling/*/SentieonDNAseq/*.vcf.gz")
//           .flatten().map{vcf -> ['SentieonDNAseq', vcf.minus(vcf.fileName)[-2].toString(), vcf]},
//         Channel.fromPath("${params.outdir}/VariantCalling/*/SentieonDNAscope/*.vcf.gz")
//           .flatten().map{vcf -> ['SentieonDNAscope', vcf.minus(vcf.fileName)[-2].toString(), vcf]},
//         Channel.fromPath("${params.outdir}/VariantCalling/*/SentieonTNscope/*.vcf.gz")
//           .flatten().map{vcf -> ['SentieonTNscope', vcf.minus(vcf.fileName)[-2].toString(), vcf]},
//         Channel.fromPath("${params.outdir}/VariantCalling/*/Strelka/*{somatic,variant}*.vcf.gz")
//           .flatten().map{vcf -> ['Strelka', vcf.minus(vcf.fileName)[-2].toString(), vcf]},
//         Channel.fromPath("${params.outdir}/VariantCalling/*/TIDDIT/*.vcf.gz")
//           .flatten().map{vcf -> ['TIDDIT', vcf.minus(vcf.fileName)[-2].toString(), vcf]}
//       ).choice(vcfToAnnotate, vcfNoAnnotate) {
//         annotate_tools == [] || (annotate_tools != [] && it[0] in annotate_tools) ? 0 : 1
//       }
//     } else if (annotate_tools == []) {
//     // Annotate user-submitted VCFs
//     // If user-submitted, Sarek assume that the idSample should be assumed automatically
//       vcfToAnnotate = Channel.fromPath(tsv_path)
//         .map{vcf -> ['userspecified', vcf.minus(vcf.fileName)[-2].toString(), vcf]}
//     } else exit 1, "specify only tools or files to annotate, not both"

//     vcfNoAnnotate.close()
//     vcfAnnotation = vcfAnnotation.mix(vcfToAnnotate)
// }

// // as now have the list of VCFs to annotate, the first step is to annotate with allele frequencies, if there are any

// (vcfSnpeff, vcfVep) = vcfAnnotation.into(2)

// vcfVep = vcfVep.map {
//   variantCaller, idSample, vcf ->
//   [variantCaller, idSample, vcf, null]
// }

// // STEP SNPEFF

// process Snpeff {
//     tag "${idSample} - ${variantCaller} - ${vcf}"

//     publishDir params.outdir, mode: params.publish_dir_mode, saveAs: {
//         if (it == "${reducedVCF}_snpEff.ann.vcf") null
//         else "Reports/${idSample}/snpEff/${it}"
//     }

//     input:
//         set variantCaller, idSample, file(vcf) from vcfSnpeff
//         file(dataDir) from ch_snpeff_cache
//         val snpeffDb from ch_snpeff_db

//     output:
//         set file("${reducedVCF}_snpEff.genes.txt"), file("${reducedVCF}_snpEff.html"), file("${reducedVCF}_snpEff.csv") into snpeffReport
//         set variantCaller, idSample, file("${reducedVCF}_snpEff.ann.vcf") into snpeffVCF

//     when: 'snpeff' in tools || 'merge' in tools

//     script:
//     reducedVCF = reduceVCF(vcf.fileName)
//     cache = (params.snpeff_cache && params.annotation_cache) ? "-dataDir \${PWD}/${dataDir}" : ""
//     """
//     snpEff -Xmx${task.memory.toGiga()}g \
//         ${snpeffDb} \
//         -csvStats ${reducedVCF}_snpEff.csv \
//         -nodownload \
//         ${cache} \
//         -canon \
//         -v \
//         ${vcf} \
//         > ${reducedVCF}_snpEff.ann.vcf

//     mv snpEff_summary.html ${reducedVCF}_snpEff.html
//     """
// }

// snpeffReport = snpeffReport.dump(tag:'snpEff report')

// // STEP COMPRESS AND INDEX VCF.1 - SNPEFF

// process CompressVCFsnpEff {
//     tag "${idSample} - ${vcf}"

//     publishDir "${params.outdir}/Annotation/${idSample}/snpEff", mode: params.publish_dir_mode

//     input:
//         set variantCaller, idSample, file(vcf) from snpeffVCF

//     output:
//         set variantCaller, idSample, file("*.vcf.gz"), file("*.vcf.gz.tbi") into (compressVCFsnpEffOut)

//     script:
//     """
//     bgzip < ${vcf} > ${vcf}.gz
//     tabix ${vcf}.gz
//     """
// }

// compressVCFsnpEffOut = compressVCFsnpEffOut.dump(tag:'VCF')

// // STEP VEP.1

// process VEP {
//     label 'VEP'
//     label 'cpus_4'

//     tag "${idSample} - ${variantCaller} - ${vcf}"

//     publishDir params.outdir, mode: params.publish_dir_mode, saveAs: {
//         if (it == "${reducedVCF}_VEP.summary.html") "Reports/${idSample}/VEP/${it}"
//         else null
//     }

//     input:
//         set variantCaller, idSample, file(vcf), file(idx) from vcfVep
//         file(dataDir) from ch_vep_cache
//         val cache_version from ch_vep_cache_version
//         file(cadd_InDels) from ch_cadd_indels
//         file(cadd_InDels_tbi) from ch_cadd_indels_tbi
//         file(cadd_WG_SNVs) from ch_cadd_wg_snvs
//         file(cadd_WG_SNVs_tbi) from ch_cadd_wg_snvs_tbi

//     output:
//         set variantCaller, idSample, file("${reducedVCF}_VEP.ann.vcf") into vepVCF
//         file("${reducedVCF}_VEP.summary.html") into vepReport

//     when: 'vep' in tools

//     script:
//     reducedVCF = reduceVCF(vcf.fileName)
//     genome = params.genome == 'smallGRCh37' ? 'GRCh37' : params.genome

//     dir_cache = (params.vep_cache && params.annotation_cache) ? " \${PWD}/${dataDir}" : "/.vep"
//     cadd = (params.cadd_cache && params.cadd_wg_snvs && params.cadd_indels) ? "--plugin CADD,whole_genome_SNVs.tsv.gz,InDels.tsv.gz" : ""
//     genesplicer = params.genesplicer ? "--plugin GeneSplicer,/opt/conda/envs/nf-core-sarek-${workflow.manifest.version}/bin/genesplicer,/opt/conda/envs/nf-core-sarek-${workflow.manifest.version}/share/genesplicer-1.0-1/human,context=200,tmpdir=\$PWD/${reducedVCF}" : "--offline"
//     """
//     mkdir ${reducedVCF}

//     vep \
//         -i ${vcf} \
//         -o ${reducedVCF}_VEP.ann.vcf \
//         --assembly ${genome} \
//         --species ${params.species} \
//         ${cadd} \
//         ${genesplicer} \
//         --cache \
//         --cache_version ${cache_version} \
//         --dir_cache ${dir_cache} \
//         --everything \
//         --filter_common \
//         --fork ${task.cpus} \
//         --format vcf \
//         --per_gene \
//         --stats_file ${reducedVCF}_VEP.summary.html \
//         --total_length \
//         --vcf

//     rm -rf ${reducedVCF}
//     """
// }

// vepReport = vepReport.dump(tag:'VEP')

// // STEP VEP.2 - VEP AFTER SNPEFF

// process VEPmerge {
//     label 'VEP'
//     label 'cpus_4'

//     tag "${idSample} - ${variantCaller} - ${vcf}"

//     publishDir params.outdir, mode: params.publish_dir_mode, saveAs: {
//         if (it == "${reducedVCF}_VEP.summary.html") "Reports/${idSample}/VEP/${it}"
//         else null
//     }

//     input:
//         set variantCaller, idSample, file(vcf), file(idx) from compressVCFsnpEffOut
//         file(dataDir) from ch_vep_cache
//         val cache_version from ch_vep_cache_version
//         file(cadd_InDels) from ch_cadd_indels
//         file(cadd_InDels_tbi) from ch_cadd_indels_tbi
//         file(cadd_WG_SNVs) from ch_cadd_wg_snvs
//         file(cadd_WG_SNVs_tbi) from ch_cadd_wg_snvs_tbi

//     output:
//         set variantCaller, idSample, file("${reducedVCF}_VEP.ann.vcf") into vepVCFmerge
//         file("${reducedVCF}_VEP.summary.html") into vepReportMerge

//     when: 'merge' in tools

//     script:
//     reducedVCF = reduceVCF(vcf.fileName)
//     genome = params.genome == 'smallGRCh37' ? 'GRCh37' : params.genome
//     dir_cache = (params.vep_cache && params.annotation_cache) ? " \${PWD}/${dataDir}" : "/.vep"
//     cadd = (params.cadd_cache && params.cadd_wg_snvs && params.cadd_indels) ? "--plugin CADD,whole_genome_SNVs.tsv.gz,InDels.tsv.gz" : ""
//     genesplicer = params.genesplicer ? "--plugin GeneSplicer,/opt/conda/envs/nf-core-sarek-${workflow.manifest.version}/bin/genesplicer,/opt/conda/envs/nf-core-sarek-${workflow.manifest.version}/share/genesplicer-1.0-1/human,context=200,tmpdir=\$PWD/${reducedVCF}" : "--offline"
//     """
//     mkdir ${reducedVCF}

//     vep \
//         -i ${vcf} \
//         -o ${reducedVCF}_VEP.ann.vcf \
//         --assembly ${genome} \
//         --species ${params.species} \
//         ${cadd} \
//         ${genesplicer} \
//         --cache \
//         --cache_version ${cache_version} \
//         --dir_cache ${dir_cache} \
//         --everything \
//         --filter_common \
//         --fork ${task.cpus} \
//         --format vcf \
//         --per_gene \
//         --stats_file ${reducedVCF}_VEP.summary.html \
//         --total_length \
//         --vcf

//     rm -rf ${reducedVCF}
//     """
// }

// vepReportMerge = vepReportMerge.dump(tag:'VEP')

// vcfCompressVCFvep = vepVCF.mix(vepVCFmerge)

// // STEP COMPRESS AND INDEX VCF.2 - VEP

// process CompressVCFvep {
//     tag "${idSample} - ${vcf}"

//     publishDir "${params.outdir}/Annotation/${idSample}/VEP", mode: params.publish_dir_mode

//     input:
//         set variantCaller, idSample, file(vcf) from vcfCompressVCFvep

//     output:
//         set variantCaller, idSample, file("*.vcf.gz"), file("*.vcf.gz.tbi") into compressVCFOutVEP

//     script:
//     """
//     bgzip < ${vcf} > ${vcf}.gz
//     tabix ${vcf}.gz
//     """
// }

// compressVCFOutVEP = compressVCFOutVEP.dump(tag:'VCF')

// /*
// ================================================================================
//                                      MultiQC
// ================================================================================
// */

// // STEP MULTIQC

// process MultiQC {
//     publishDir "${params.outdir}/Reports/MultiQC", mode: params.publish_dir_mode

//     input:
//         file (multiqcConfig) from multiqc_config
//         file (mqc_custom_config) from multiqc_custom_config.collect().ifEmpty([])
//         file (versions) from ch_software_versions_yaml.collect()
//         file workflow_summary from workflow_summary.collectFile(name: "workflow_summary_mqc.yaml")
//         file ('bamQC/*') from bamQCReport.collect().ifEmpty([])
//         file ('BCFToolsStats/*') from bcftoolsReport.collect().ifEmpty([])
//         file ('FastQC/*') from fastQCReport.collect().ifEmpty([])
//         file ('TrimmedFastQC/*') from trimGaloreReport.collect().ifEmpty([])
//         file ('MarkDuplicates/*') from duplicates_marked_report.collect().ifEmpty([])
//         file ('DuplicatesMarked/*.recal.table') from baseRecalibratorReport.collect().ifEmpty([])
//         file ('SamToolsStats/*') from samtoolsStatsReport.collect().ifEmpty([])
//         file ('snpEff/*') from snpeffReport.collect().ifEmpty([])
//         file ('VCFTools/*') from vcftoolsReport.collect().ifEmpty([])

//     output:
//         file "*multiqc_report.html" into ch_multiqc_report
//         file "*_data"
//         file "multiqc_plots"

//     when: !('multiqc' in skipQC)

//     script:
//     rtitle = custom_runName ? "--title \"$custom_runName\"" : ''
//     rfilename = custom_runName ? "--filename " + custom_runName.replaceAll('\\W','_').replaceAll('_+','_') + "_multiqc_report" : ''
//     custom_config_file = params.multiqc_config ? "--config $mqc_custom_config" : ''
//     """
//     multiqc -f ${rtitle} ${rfilename} ${custom_config_file} .
//     """
// }

// ch_multiqc_report.dump(tag:'MultiQC')
