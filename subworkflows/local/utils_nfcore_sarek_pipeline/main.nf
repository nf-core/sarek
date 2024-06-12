//
// Subworkflow with functionality specific to the nf-core/sarek pipeline
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { paramsSummaryMap          } from 'plugin/nf-validation'
include { fromSamplesheet           } from 'plugin/nf-validation'
include { UTILS_NEXTFLOW_PIPELINE   } from '../../nf-core/utils_nextflow_pipeline'
include { UTILS_NFVALIDATION_PLUGIN } from '../../nf-core/utils_nfvalidation_plugin'
include { UTILS_NFCORE_PIPELINE     } from '../../nf-core/utils_nfcore_pipeline'
include { completionEmail           } from '../../nf-core/utils_nfcore_pipeline'
include { completionSummary         } from '../../nf-core/utils_nfcore_pipeline'
include { dashedLine                } from '../../nf-core/utils_nfcore_pipeline'
include { getWorkflowVersion        } from '../../nf-core/utils_nfcore_pipeline'
include { imNotification            } from '../../nf-core/utils_nfcore_pipeline'
include { logColours                } from '../../nf-core/utils_nfcore_pipeline'
include { MULTIQC                   } from '../../../modules/nf-core/multiqc'
include { paramsSummaryMultiqc      } from '../../nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML    } from '../../nf-core/utils_nfcore_pipeline'
include { workflowCitation          } from '../../nf-core/utils_nfcore_pipeline'
include { SAMPLESHEET_TO_CHANNEL    } from '../samplesheet_to_channel'

/*
========================================================================================
    SUBWORKFLOW TO INITIALISE PIPELINE
========================================================================================
*/

workflow PIPELINE_INITIALISATION {

    take:
    version           // boolean: Display version and exit
    help              // boolean: Display help text
    validate_params   // boolean: Boolean whether to validate parameters against the schema at runtime
    monochrome_logs   // boolean: Do not use coloured log outputs
    nextflow_cli_args //   array: List of positional nextflow CLI args
    outdir            //  string: The output directory where the results will be saved
    input             //  string: Path to input samplesheet

    main:

    versions = Channel.empty()

    //
    // Print version and exit if required and dump pipeline parameters to JSON file
    //
    UTILS_NEXTFLOW_PIPELINE (
        version,
        true,
        outdir,
        workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1
    )

    //
    // Validate parameters and generate parameter summary to stdout
    //
    pre_help_text = nfCoreLogo(monochrome_logs)
    post_help_text = '\n' + workflowCitation() + '\n' + dashedLine(monochrome_logs)
    def String workflow_command = "nextflow run ${workflow.manifest.name} -profile <docker/singularity/.../institute> --input samplesheet.csv --outdir <OUTDIR>"
    UTILS_NFVALIDATION_PLUGIN (
        help,
        workflow_command,
        pre_help_text,
        post_help_text,
        validate_params,
        "nextflow_schema.json"
    )

    //
    // Check config provided to the pipeline
    //
    UTILS_NFCORE_PIPELINE(nextflow_cli_args)

    //
    // Custom validation for pipeline parameters
    //
    validateInputParameters()

    // Check input path parameters to see if they exist
    def checkPathParamList = [
        params.ascat_alleles,
        params.ascat_loci,
        params.ascat_loci_gc,
        params.ascat_loci_rt,
        params.bwa,
        params.bwamem2,
        params.bcftools_annotations,
        params.bcftools_annotations_tbi,
        params.bcftools_header_lines,
        params.cf_chrom_len,
        params.chr_dir,
        params.cnvkit_reference,
        params.dbnsfp,
        params.dbnsfp_tbi,
        params.dbsnp,
        params.dbsnp_tbi,
        params.dict,
        params.dragmap,
        params.fasta,
        params.fasta_fai,
        params.germline_resource,
        params.germline_resource_tbi,
        params.input,
        params.intervals,
        params.known_indels,
        params.known_indels_tbi,
        params.known_snps,
        params.known_snps_tbi,
        params.mappability,
        params.multiqc_config,
        params.ngscheckmate_bed,
        params.pon,
        params.pon_tbi,
        params.sentieon_dnascope_model,
        params.spliceai_indel,
        params.spliceai_indel_tbi,
        params.spliceai_snv,
        params.spliceai_snv_tbi
    ]

// only check if we are using the tools
if (params.tools && (params.tools.split(',').contains('snpeff') || params.tools.split(',').contains('merge'))) checkPathParamList.add(params.snpeff_cache)
if (params.tools && (params.tools.split(',').contains('vep')    || params.tools.split(',').contains('merge'))) checkPathParamList.add(params.vep_cache)

    // def retrieveInput(need_input, step, outdir) {

    params.input_restart = retrieveInput((!params.build_only_index && !params.input), params.step, params.outdir)

    ch_from_samplesheet = params.build_only_index ? Channel.empty() : params.input ? Channel.fromSamplesheet("input") : Channel.fromSamplesheet("input_restart")

    SAMPLESHEET_TO_CHANNEL(
        ch_from_samplesheet,
        params.aligner,
        params.ascat_alleles,
        params.ascat_loci,
        params.ascat_loci_gc,
        params.ascat_loci_rt,
        params.bcftools_annotations,
        params.bcftools_annotations_tbi,
        params.bcftools_header_lines,
        params.build_only_index,
        params.dbsnp,
        params.fasta,
        params.germline_resource,
        params.intervals,
        params.joint_germline,
        params.joint_mutect2,
        params.known_indels,
        params.known_snps,
        params.no_intervals,
        params.pon,
        params.sentieon_dnascope_emit_mode,
        params.sentieon_haplotyper_emit_mode,
        params.seq_center,
        params.seq_platform,
        params.skip_tools,
        params.step,
        params.tools,
        params.umi_read_structure,
        params.wes)

    emit:
    samplesheet = SAMPLESHEET_TO_CHANNEL.out.input_sample
    versions
}

/*
========================================================================================
    SUBWORKFLOW FOR VERSIONS AND REPORT GATHERING AND MULTIQC
========================================================================================
*/

workflow GATHER_REPORTS_VERSIONS {
    take:
    outdir
    multiqc_config
    multiqc_logo
    multiqc_methods_description
    versions
    reports

    main:

    //
    // Collate and save software versions
    //
    version_yaml = Channel.empty()
    version_yaml = softwareVersionsToYAML(versions)
        .collectFile(storeDir: "${outdir}/pipeline_info", name: 'nf_core_sarek_software_mqc_versions.yml', sort: true, newLine: true)

    //
    // MODULE: MultiQC
    //
    ch_multiqc_files                      = Channel.empty()
    ch_multiqc_config                     = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config              = multiqc_config ? Channel.fromPath(multiqc_config, checkIfExists: true) : Channel.empty()
    ch_multiqc_logo                       = multiqc_logo ? Channel.fromPath(multiqc_logo, checkIfExists: true) : Channel.empty()
    summary_params                        = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary                   = Channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_custom_methods_description = multiqc_methods_description ? file(multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = Channel.value(methodsDescriptionText(ch_multiqc_custom_methods_description))
    ch_multiqc_files                      = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files                      = ch_multiqc_files.mix(version_yaml)
    ch_multiqc_files                      = ch_multiqc_files.mix(reports)
    ch_multiqc_files                      = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml', sort: true))

    MULTIQC(
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList()
    )

    emit:
    report = MULTIQC.out.report.toList()

}

/*
========================================================================================
    SUBWORKFLOW FOR PIPELINE COMPLETION
========================================================================================
*/

workflow PIPELINE_COMPLETION {

    take:
    email           //  string: email address
    email_on_fail   //  string: email address sent on pipeline failure
    plaintext_email // boolean: Send plain-text email instead of HTML
    outdir          //    path: Path to output directory where results will be published
    monochrome_logs // boolean: Disable ANSI colour codes in log output
    hook_url        //  string: hook URL for notifications
    multiqc_report  //  string: Path to MultiQC report

    main:

    summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")

    //
    // Completion email and summary
    //
    workflow.onComplete {
        if (email || email_on_fail) {
            completionEmail(summary_params, email, email_on_fail, plaintext_email, outdir, monochrome_logs, multiqc_report.toList())
        }

        completionSummary(monochrome_logs)

        if (hook_url) {
            imNotification(summary_params, hook_url)
        }
    }

    workflow.onError {
        log.error "Pipeline failed. Please refer to troubleshooting docs: https://nf-co.re/docs/usage/troubleshooting"
    }
}

/*
========================================================================================
    FUNCTIONS
========================================================================================
*/
//
// Check and validate pipeline parameters
//
def validateInputParameters() {
    genomeExistsError()
}

//
// Validate channels from input samplesheet
//
def validateInputSamplesheet(input) {
    def (metas, fastqs) = input[1..2]

    // Check that multiple runs of the same sample are of the same datatype i.e. single-end / paired-end
    def endedness_ok = metas.collect{ it.single_end }.unique().size == 1
    if (!endedness_ok) {
        error("Please check input samplesheet -> Multiple runs of a sample must be of the same datatype i.e. single-end or paired-end: ${metas[0].id}")
    }

    return [ metas[0], fastqs ]
}

//
// Exit pipeline if incorrect --genome key provided
//
def genomeExistsError() {
    if (params.genomes && params.genome && !params.genomes.containsKey(params.genome)) {
        def error_string = "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
            "  Genome '${params.genome}' not found in any config files provided to the pipeline.\n" +
            "  Currently, the available genome keys are:\n" +
            "  ${params.genomes.keySet().join(", ")}\n" +
            "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
        error(error_string)
    }
}

//
// Generate methods description for MultiQC
//
def toolCitationText() {
    // TODO nf-core: Optionally add in-text citation tools to this list.
    // Can use ternary operators to dynamically construct based conditions, e.g. params["run_xyz"] ? "Tool (Foo et al. 2023)" : "",
    // Uncomment function in methodsDescriptionText to render in MultiQC report
    def citation_text = [
            "Tools used in the workflow included:",
            "FastQC (Andrews 2010),",
            "MultiQC (Ewels et al. 2016)",
            "."
        ].join(' ').trim()

    return citation_text
}

def toolBibliographyText() {
    // TODO nf-core: Optionally add bibliographic entries to this list.
    // Can use ternary operators to dynamically construct based conditions, e.g. params["run_xyz"] ? "<li>Author (2023) Pub name, Journal, DOI</li>" : "",
    // Uncomment function in methodsDescriptionText to render in MultiQC report
    def reference_text = [
            "<li>Andrews S, (2010) FastQC, URL: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/).</li>",
            "<li>Ewels, P., Magnusson, M., Lundin, S., & Käller, M. (2016). MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics , 32(19), 3047–3048. doi: /10.1093/bioinformatics/btw354</li>"
        ].join(' ').trim()

    return reference_text
}

def methodsDescriptionText(mqc_methods_yaml) {
    // Convert  to a named map so can be used as with familar NXF ${workflow} variable syntax in the MultiQC YML file
    def meta = [:]
    meta.workflow = workflow.toMap()
    meta["manifest_map"] = workflow.manifest.toMap()

    // Pipeline DOI
    if (meta.manifest_map.doi) {
        // Using a loop to handle multiple DOIs
        // Removing `https://doi.org/` to handle pipelines using DOIs vs DOI resolvers
        // Removing ` ` since the manifest.doi is a string and not a proper list
        def temp_doi_ref = ""
        String[] manifest_doi = meta.manifest_map.doi.tokenize(",")
        for (String doi_ref: manifest_doi) temp_doi_ref += "(doi: <a href=\'https://doi.org/${doi_ref.replace("https://doi.org/", "").replace(" ", "")}\'>${doi_ref.replace("https://doi.org/", "").replace(" ", "")}</a>), "
        meta["doi_text"] = temp_doi_ref.substring(0, temp_doi_ref.length() - 2)
    } else meta["doi_text"] = ""
    meta["nodoi_text"] = meta.manifest_map.doi ? "" : "<li>If available, make sure to update the text to include the Zenodo DOI of the pipeline version used. </li>"

    // Tool references
    meta["tool_citations"] = ""
    meta["tool_bibliography"] = ""

    // TODO nf-core: Only uncomment below if logic in toolCitationText/toolBibliographyText has been filled!
    // meta["tool_citations"] = toolCitationText().replaceAll(", \\.", ".").replaceAll("\\. \\.", ".").replaceAll(", \\.", ".")
    // meta["tool_bibliography"] = toolBibliographyText()


    def methods_text = mqc_methods_yaml.text

    def engine =  new groovy.text.SimpleTemplateEngine()
    def description_html = engine.createTemplate(methods_text).make(meta)

    return description_html.toString()
}

//
// nf-core/sarek logo
//
def nfCoreLogo(monochrome_logs=true) {
    Map colors = logColours(monochrome_logs)
    String.format(
        """\n
        ${dashedLine(monochrome_logs)}
                                                ${colors.green},--.${colors.black}/${colors.green},-.${colors.reset}
        ${colors.blue}        ___     __   __   __   ___     ${colors.green}/,-._.--~\'${colors.reset}
        ${colors.blue}  |\\ | |__  __ /  ` /  \\ |__) |__         ${colors.yellow}}  {${colors.reset}
        ${colors.blue}  | \\| |       \\__, \\__/ |  \\ |___     ${colors.green}\\`-._,-`-,${colors.reset}
                                                ${colors.green}`._,._,\'${colors.reset}
        ${colors.white}      ____${colors.reset}
        ${colors.white}    .´ _  `.${colors.reset}
        ${colors.white}   /  ${colors.green}|\\${colors.reset}`-_ \\${colors.reset}     ${colors.blue} __        __   ___     ${colors.reset}
        ${colors.white}  |   ${colors.green}| \\${colors.reset}  `-|${colors.reset}    ${colors.blue}|__`  /\\  |__) |__  |__/${colors.reset}
        ${colors.white}   \\ ${colors.green}|   \\${colors.reset}  /${colors.reset}     ${colors.blue}.__| /¯¯\\ |  \\ |___ |  \\${colors.reset}
        ${colors.white}    `${colors.green}|${colors.reset}____${colors.green}\\${colors.reset}´${colors.reset}

        ${colors.purple}  ${workflow.manifest.name} ${getWorkflowVersion()}${colors.reset}
        ${dashedLine(monochrome_logs)}
        """.stripIndent()
    )
}

//
// retrieveInput
//
def retrieveInput(need_input, step, outdir) {
    def input = null
    if (!params.input && !params.build_only_index) {
        switch (step) {
            case 'mapping':                 Nextflow.error("Can't start with step $step without samplesheet")
                                            break
            case 'markduplicates':          log.warn("Using file ${outdir}/csv/mapped.csv");
                                            input = outdir + "/csv/mapped.csv"
                                            break
            case 'prepare_recalibration':   log.warn("Using file ${outdir}/csv/markduplicates_no_table.csv");
                                            input = outdir + "/csv/markduplicates_no_table.csv"
                                            break
            case 'recalibrate':             log.warn("Using file ${outdir}/csv/markduplicates.csv");
                                            input = outdir + "/csv/markduplicates.csv"
                                            break
            case 'variant_calling':         log.warn("Using file ${outdir}/csv/recalibrated.csv");
                                            input = outdir + "/csv/recalibrated.csv"
                                            break
            // case 'controlfreec':         csv_file = file("${outdir}/variant_calling/csv/control-freec_mpileup.csv", checkIfExists: true); break
            case 'annotate':                log.warn("Using file ${outdir}/csv/variantcalled.csv");
                                            input = outdir + "/csv/variantcalled.csv"
                                            break
            default:                        log.warn("Please provide an input samplesheet to the pipeline e.g. '--input samplesheet.csv'")
                                            Nextflow.error("Unknown step $step")
        }
    }
    return input
}
