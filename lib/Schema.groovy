/*
 * This file holds several functions used to perform JSON parameter validation, help and summary rendering for the nf-core pipeline template.
 */

import groovy.json.JsonSlurper

class JSON {
    /*
     * This method tries to read a JSON params file
     */
    private static LinkedHashMap params_get(String path) {
        def usage = new LinkedHashMap()
        try {
            usage = params_try(path)
        } catch (Exception e) {
            println "Could not read parameters settings from JSON. $e"
            usage = new LinkedHashMap()
        }
        return usage
    }

    /*
    Method to actually read in JSON file using Groovy.
    Group (as Key), values are all parameters
        - Parameter1 as Key, Description as Value
        - Parameter2 as Key, Description as Value
        ....
    Group
        -
    */
    private static LinkedHashMap params_try(String path) throws Exception {

        def json = new File(path).text
        def Map usage = (Map) new JsonSlurper().parseText(json).get('properties')

        /* Tree looks like this in nf-core schema
        *  properties <- this is what the first get('properties') gets us
             group 1
               properties
               description
             group 2
               properties
               description
             group 3
               properties
               description
        */
        def output_map = new LinkedHashMap()

        // Lets go deeper
        usage.each { key, val ->
            def Map submap = usage."$key".properties // Gets the property object of the group
            def sub_params = new LinkedHashMap()
            submap.each { innerkey, value ->
                sub_params.put("$innerkey", "$value.description")
            }
            output_map.put("$key", sub_params)
        }
        return output_map
    }

    static String params_help(path, command) {
          String output = "Typical pipeline command:\n\n"
          output += "    ${command}\n\n"
          output += params_beautify(params_get(path))
    }

    static String params_beautify(usage) {
        String output = ""
        for (group in usage.keySet()) {
            output += group + "\n"
            def params = usage.get(group)  // This gets the parameters of that particular group
            for (par in params.keySet()) {
                output+= "    \u001B[1m" + par.padRight(27) + "\u001B[1m" + params.get(par) + "\n"
            }
            output += "\n"
        }
        return output
    }

    private static LinkedHashMap params_summary(workflow, params, run_name) {
        def Map summary = [:]
        if (workflow.revision) summary['Pipeline Release'] = workflow.revision
        summary['Run Name']         = run_name ?: workflow.runName
        // TODO nf-core: Report custom parameters here
        summary['Input']            = params.input
        summary['Fasta File']       = params.fasta
        summary['Max Resources']    = "$params.max_memory memory, $params.max_cpus cpus, $params.max_time time per job"
        if (workflow.containerEngine) summary['Container'] = "$workflow.containerEngine - $workflow.container"
        summary['Output dir']       = params.outdir
        summary['Launch dir']       = workflow.launchDir
        summary['Working dir']      = workflow.workDir
        summary['Script dir']       = workflow.projectDir
        summary['User']             = workflow.userName
        if (workflow.profile.contains('awsbatch')) {
            summary['AWS Region']   = params.awsregion
            summary['AWS Queue']    = params.awsqueue
            summary['AWS CLI']      = params.awscli
        }
        summary['Config Profile'] = workflow.profile
        if (params.config_profile_description) summary['Config Profile Descr']   = params.config_profile_description
        if (params.config_profile_contact)     summary['Config Profile Contact'] = params.config_profile_contact
        if (params.config_profile_url)         summary['Config Profile URL']     = params.config_profile_url
        summary['Config Files'] = workflow.configFiles.join(', ')
        if (params.email || params.email_on_fail) {
            summary['E-mail Address']    = params.email
            summary['E-mail on failure'] = params.email_on_fail
            summary['MultiQC maxsize']   = params.max_multiqc_email_size
        }
        return summary
    }

    static String params_mqc_summary(summary) {
        String yaml_file_text  = """
        id: 'nf-core-sarek-summary'
        description: " - this information is collected when the pipeline is started."
        section_name: 'nf-core/sarek Workflow Summary'
        section_href: 'https://github.com/nf-core/sarek'
        plot_type: 'html'
        data: |
            <dl class=\"dl-horizontal\">
            ${summary.collect { k,v -> "            <dt>$k</dt><dd><samp>${v ?: '<span style=\"color:#999999;\">N/A</a>'}</samp></dd>" }.join("\n")}
            </dl>
        """.stripIndent()

        return yaml_file_text
    }
}
