/*
 * Output Markdown documentation to HTML
 */
process OUTPUT_DOCUMENTATION {
    publishDir "${params.outdir}/pipeline_info", mode: params.publish_dir_mode

    container "nfcore/sarek:dsl2"

    conda (params.conda ? "$baseDir/environment.yml" : null)

    input:
        path output_docs
        path images

    output:
        path "results_description.html"

    script:
    """
    markdown_to_html.py ${output_docs} -o results_description.html
    """
}
