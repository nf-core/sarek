//
// Subworkflow with functionality specific to the nf-core/sarek pipeline
//

/*
========================================================================================
    IMPORT MODULES/SUBWORKFLOWS
========================================================================================
*/

include { UTILS_NEXTFLOW_PIPELINE    } from '../../nf-core/utils_nextflow_pipeline'
include { UTILS_NFCORE_PIPELINE      } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NFVALIDATION_PLUGIN  } from '../../nf-core/utils_nfvalidation_plugin'
include { completionEmail            } from '../../nf-core/utils_nfcore_pipeline'
include { completionSummary          } from '../../nf-core/utils_nfcore_pipeline'
include { dashedLine                 } from '../../nf-core/utils_nfcore_pipeline'
include { getWorkflowVersion         } from '../../nf-core/utils_nextflow_pipeline'
include { imNotification             } from '../../nf-core/utils_nfcore_pipeline'
include { logColours                 } from '../../nf-core/utils_nfcore_pipeline'
include { workflowCitation           } from '../../nf-core/utils_nfcore_pipeline'

/*
========================================================================================
    SUBWORKFLOW TO INITIALISE PIPELINE
========================================================================================
*/

workflow PIPELINE_INITIALISATION {

    main:

    //
    // Print version and exit if required and dump pipeline parameters to JSON file
    //
    UTILS_NEXTFLOW_PIPELINE (
        params.version,
        true,
        params.outdir,
        workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1
    )

    //
    // Validate parameters and generate parameter summary to stdout
    //
    def pre_help_text = nfCoreLogo(getWorkflowVersion())
    def post_help_text = '\n' + workflowCitation() + '\n' + dashedLine()
    def String workflow_command = "nextflow run ${workflow.manifest.name} -profile <docker/singularity/.../institute> --input samplesheet.csv --outdir <OUTDIR>"
    UTILS_NFVALIDATION_PLUGIN (
        params.help,
        workflow_command,
        pre_help_text,
        post_help_text,
        params.validate_params,
        "nextflow_schema.json"
    )

    //
    // Check config provided to the pipeline
    //
    UTILS_NFCORE_PIPELINE ()

    emit:
    summary_params = UTILS_NFVALIDATION_PLUGIN.out.summary_params
}

/*
========================================================================================
    SUBWORKFLOW FOR PIPELINE COMPLETION
========================================================================================
*/

workflow PIPELINE_COMPLETION {

    take:
    versions       // channel: software tools versions
    email          //  string: email address
    email_on_fail  //  string: email address sent on pipeline failure
    hook_url       //  string: hook URL for notifications
    summary_params //     map: Groovy map of the parameters used in the pipeline

    main:

    //
    // Completion email and summary
    //
    workflow.onComplete {
        if (email || email_on_fail) {
            completionEmail(summary_params)
        }

        completionSummary()

        if (hook_url) {
            imNotification(summary_params)
        }

    }
}

/*
========================================================================================
    FUNCTIONS
========================================================================================
*/

def getGenomeAttribute(params, attribute) {
    if (params.genomes && params.genome && params.genomes.containsKey(params.genome)) {
        if (params.genomes[ params.genome ].containsKey(attribute)) {
            return params.genomes[ params.genome ][ attribute ]
        }
    }
    return null
}

def nfCoreLogo(workflow_version) {
    Map colors = logColours()
    String.format(
        """\n
        ${dashedLine()}
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
        ${colors.purple}  ${workflow.manifest.name} ${workflow_version}${colors.reset}
        ${dashedLine()}
        """.stripIndent()
    )
}


def retrieveInput(params, log){
    def input = null
    if (!params.input && !params.build_only_index) {
        switch (params.step) {
            case 'mapping':                 Nextflow.error("Can't start with step $params.step without samplesheet")
                                            break
            case 'markduplicates':          log.warn("Using file ${params.outdir}/csv/mapped.csv");
                                            input = params.outdir + "/csv/mapped.csv"
                                            break
            case 'prepare_recalibration':   log.warn("Using file ${params.outdir}/csv/markduplicates_no_table.csv");
                                            input = params.outdir + "/csv/markduplicates_no_table.csv"
                                            break
            case 'recalibrate':             log.warn("Using file ${params.outdir}/csv/markduplicates.csv");
                                            input = params.outdir + "/csv/markduplicates.csv"
                                            break
            case 'variant_calling':         log.warn("Using file ${params.outdir}/csv/recalibrated.csv");
                                            input = params.outdir + "/csv/recalibrated.csv"
                                            break
            // case 'controlfreec':         csv_file = file("${params.outdir}/variant_calling/csv/control-freec_mpileup.csv", checkIfExists: true); break
            case 'annotate':                log.warn("Using file ${params.outdir}/csv/variantcalled.csv");
                                            input = params.outdir + "/csv/variantcalled.csv"
                                            break
            default:                        log.warn("Please provide an input samplesheet to the pipeline e.g. '--input samplesheet.csv'")
                                            Nextflow.error("Unknown step $params.step")
        }
    }
    return input
}
