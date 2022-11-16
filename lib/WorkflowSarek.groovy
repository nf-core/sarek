//
// This file holds several functions specific to the workflow/sarek.nf in the nf-core/sarek pipeline
//

import groovy.text.SimpleTemplateEngine

class WorkflowSarek {

    //
    // Check and validate parameters
    //
    public static void initialise(params, log) {
        genomeExistsError(params, log)


        if (!params.fasta) {
            log.error "Genome fasta file not specified with e.g. '--fasta genome.fa' or via a detectable config file."
            System.exit(1)
        }
    }

    //
    // Get workflow summary for MultiQC
    //
    public static String paramsSummaryMultiqc(workflow, summary) {
        String summary_section = ''
        for (group in summary.keySet()) {
            def group_params = summary.get(group)  // This gets the parameters of that particular group
            if (group_params) {
                summary_section += "    <p style=\"font-size:110%\"><b>$group</b></p>\n"
                summary_section += "    <dl class=\"dl-horizontal\">\n"
                for (param in group_params.keySet()) {
                    summary_section += "        <dt>$param</dt><dd><samp>${group_params.get(param) ?: '<span style=\"color:#999999;\">N/A</a>'}</samp></dd>\n"
                }
                summary_section += "    </dl>\n"
            }
        }

        String yaml_file_text  = "id: '${workflow.manifest.name.replace('/','-')}-summary'\n"
        yaml_file_text        += "description: ' - this information is collected when the pipeline is started.'\n"
        yaml_file_text        += "section_name: '${workflow.manifest.name} Workflow Summary'\n"
        yaml_file_text        += "section_href: 'https://github.com/${workflow.manifest.name}'\n"
        yaml_file_text        += "plot_type: 'html'\n"
        yaml_file_text        += "data: |\n"
        yaml_file_text        += "${summary_section}"
        return yaml_file_text
    }

    public static String methodsDescriptionText(run_workflow, mqc_methods_yaml) {
        // Convert  to a named map so can be used as with familar NXF ${workflow} variable syntax in the MultiQC YML file
        def meta = [:]
        meta.workflow = run_workflow.toMap()
        meta["manifest_map"] = run_workflow.manifest.toMap()

        meta["doi_text"] = meta.manifest_map.doi ? "(doi: <a href=\'https://doi.org/${meta.manifest_map.doi}\'>${meta.manifest_map.doi}</a>)" : ""
        meta["nodoi_text"] = meta.manifest_map.doi ? "": "<li>If available, make sure to update the text to include the Zenodo DOI of version of the pipeline used. </li>"

        def methods_text = mqc_methods_yaml.text

        def engine =  new SimpleTemplateEngine()
        def description_html = engine.createTemplate(methods_text).make(meta)

        return description_html
    }//
    // Exit pipeline if incorrect --genome key provided
    //
    private static void genomeExistsError(params, log) {
        if (params.genomes && params.genome && !params.genomes.containsKey(params.genome)) {
            log.error "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
                "  Genome '${params.genome}' not found in any config files provided to the pipeline.\n" +
                "  Currently, the available genome keys are:\n" +
                "  ${params.genomes.keySet().join(", ")}\n" +
                "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
            System.exit(1)
        }
    }

    public static String retrieveInput(params, log){
        switch (params.step) {
            case 'mapping':                 log.warn "Can't start with step $params.step without samplesheet"
                                            System.exit(1);
                                            break
            case 'markduplicates':          log.warn "Using file ${params.outdir}/csv/mapped.csv"
                                            params.putIfAbsent("input","${params.outdir}/csv/mapped.csv");
                                            break
            case 'prepare_recalibration':   log.warn "Using file ${params.outdir}/csv/markduplicates_no_table.csv"
                                            params.putIfAbsent("input", "${params.outdir}/csv/markduplicates_no_table.csv");
                                            break
            case 'recalibrate':             log.warn "Using file ${params.outdir}/csv/markduplicates.csv"
                                            params.putIfAbsent("input", "${params.outdir}/csv/markduplicates.csv");
                                            break
            case 'variant_calling':         log.warn "Using file ${params.outdir}/csv/recalibrated.csv"
                                            params.putIfAbsent("input", "${params.outdir}/csv/recalibrated.csv");
                                            break
            // case 'controlfreec':         csv_file = file("${params.outdir}/variant_calling/csv/control-freec_mpileup.csv", checkIfExists: true); break
            case 'annotate':                log.warn "Using file ${params.outdir}/csv/variantcalled.csv"
                                            params.putIfAbsent("input","${params.outdir}/csv/variantcalled.csv");
                                            break
            default:    log.warn "Please provide an input samplesheet to the pipeline e.g. '--input samplesheet.csv'"
                        exit 1, "Unknown step $params.step"
        }
    }
}
