//
// This file holds several functions specific to the main.nf workflow in the nf-core/sarek pipeline
//

class WorkflowMain {

    //
    // Citation string for pipeline
    //
    public static String citation(workflow) {
        return "If you use ${workflow.manifest.name} for your analysis please cite:\n\n" +
            "* The pipeline\n" +
            "  https://doi.org/10.12688/f1000research.16665.2\n" +
            "  https://doi.org/10.5281/zenodo.4468605\n\n" +
            "* The nf-core framework\n" +
            "  https://doi.org/10.1038/s41587-020-0439-x\n\n" +
            "* Software dependencies\n" +
            "  https://github.com/${workflow.manifest.name}/blob/master/CITATIONS.md"
    }

    //
    // Print help to screen if required
    //
    public static String help(workflow, params, log) {
        def command = "nextflow run ${workflow.manifest.name} --input samplesheet.tsv --genome GRCh38 -profile docker"
        def help_string = ''
        help_string += NfcoreTemplate.logo(workflow, params.monochrome_logs)
        help_string += NfcoreSchema.paramsHelp(workflow, params, command)
        help_string += '\n' + citation(workflow) + '\n'
        help_string += NfcoreTemplate.dashedLine(params.monochrome_logs)
        return help_string
    }

    //
    // Print parameter summary log to screen
    //
    public static String paramsSummaryLog(workflow, params, log) {
        def summary_log = ''
        summary_log += NfcoreTemplate.logo(workflow, params.monochrome_logs)
        summary_log += NfcoreSchema.paramsSummaryLog(workflow, params)
        summary_log += '\n' + citation(workflow) + '\n'
        summary_log += NfcoreTemplate.dashedLine(params.monochrome_logs)
        return summary_log
    }

    //
    // Validate parameters and print summary to screen
    //
    public static void initialise(workflow, params, log) {
        // Print help to screen if required
        if (params.help) {
            log.info help(workflow, params, log)
            System.exit(0)
        }

        // Validate workflow parameters via the JSON schema
        if (params.validate_params) {
            NfcoreSchema.validateParameters(workflow, params, log)
        }

        // Print parameter summary log to screen
        log.info paramsSummaryLog(workflow, params, log)

        // Check that a -profile or Nextflow config has been provided to run the pipeline
        NfcoreTemplate.checkConfigProvided(workflow, log)

        // Check that conda channels are set-up correctly
        if (params.enable_conda) {
            Utils.checkCondaChannels(log)
        }

        // Check AWS batch settings
        NfcoreTemplate.awsBatch(workflow, params)

        // Check input has been provided
        if (!params.input) {
            log.warn "No samplesheet specified, attempting to restart from csv files present in ${params.outdir}"

            // switch (params.step) {
            //     case 'mapping': log.warn "Can't start with step $params.step without samplesheet"
            //                     System.exit(1);
            //                     break
            //     //case 'markduplicates':          log.warn "Using file ${params.outdir}/preprocessing/csv/markduplicates_no_table.csv"
            //     //                                params.input = "${params.outdir}/preprocessing/csv/markduplicates_no_table.csv";
            //     //                                break
            //     case 'prepare_recalibration':   log.warn "Using file ${params.outdir}/preprocessing/csv/markduplicates_no_table.csv"
            //                                     params.input = "${params.outdir}/preprocessing/csv/markduplicates_no_table.csv";
            //                                     break
            //     case 'recalibrate':             log.warn "Using file ${params.outdir}/preprocessing/csv/markduplicates.csv"
            //                                     params.input = "${params.outdir}/preprocessing/csv/markduplicates.csv";
            //                                     break
            //     case 'variant_calling':         log.warn "Using file ${params.outdir}/preprocessing/csv/recalibrated.csv"
            //                                     params.input = "${params.outdir}/preprocessing/csv/recalibrated.csv";
            //                                     break
            // //    // case 'controlfreec':         csv_file = file("${params.outdir}/variant_calling/csv/control-freec_mpileup.csv", checkIfExists: true); break
            // //    case 'annotate':              csv_file = file("${params.outdir}/variant_calling/csv/recalibrated.csv",          checkIfExists: true); break
            //     default:    log.warn "Please provide an input samplesheet to the pipeline e.g. '--input samplesheet.csv'"
            //                 exit 1, "Unknown step $params.step"

            //log.warn "Please provide an input samplesheet to the pipeline e.g. '--input samplesheet.csv'"
            System.exit(1)
        }
    }

    //
    // Get attribute from genome config file e.g. fasta
    //
    public static String getGenomeAttribute(params, attribute) {
        def val = ''
        if (params.genomes && params.genome && params.genomes.containsKey(params.genome)) {
            if (params.genomes[ params.genome ].containsKey(attribute)) {
                val = params.genomes[ params.genome ][ attribute ]
            }
        }
        return val
    }
}
