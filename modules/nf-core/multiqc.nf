// Has the run name been specified by the user?
// this has the bonus effect of catching both -name and --name
custom_runName = params.name
if (!(workflow.runName ==~ /[a-z]+_[a-z]+/)) {
    custom_runName = workflow.runName
}

/*
 * MultiQC
 */
process MULTIQC {
    publishDir "${params.outdir}/multiqc", mode: params.publish_dir_mode

    input:
    path multiqc_config
    path mqc_custom_config
    // TODO nf-core: Add in log files from your new processes for MultiQC to find!
    path fastqc
    path software_versions
    val workflow_summary

    output:
    path "*multiqc_report.html"
    path "*_data"
    path "multiqc_plots"

    script:
    rtitle = custom_runName ? "--title \"$custom_runName\"" : ''
    rfilename = custom_runName ? "--filename " + custom_runName.replaceAll('\\W','_').replaceAll('_+','_') + "_multiqc_report" : ''
    custom_config_file = params.multiqc_config ? "--config $mqc_custom_config" : ''
    // TODO nf-core: Specify which MultiQC modules to use with -m for a faster run time
    """
    echo '$workflow_summary' > workflow_summary_mqc.yaml
    multiqc -f $rtitle $rfilename $custom_config_file .
    """
}
