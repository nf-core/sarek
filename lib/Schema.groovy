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

    private static LinkedHashMap params_summary(workflow, params, run_name, step, tools, skip_qc, annotate_tools) {
        def Map summary = [:]
        if (workflow.revision) summary['Pipeline Release'] = workflow.revision
        summary['Run Name']         = run_name ?: workflow.runName
        summary['Max Resources']    = "${params.max_memory} memory, ${params.max_cpus} cpus, ${params.max_time} time per job"
        if (workflow.containerEngine) summary['Container'] = "${workflow.containerEngine} - ${workflow.container}"
        summary['Input']             = params.input
        summary['Step']              = step
        summary['Genome']            = params.genome
        if (params.no_intervals && step != 'annotate')  summary['Intervals']         = 'Do not use'
        summary['Nucleotides/s']     = params.nucleotides_per_second
        if (params.sentieon)            summary['Sention']                           = "Using Sentieon for Preprocessing and/or Variant Calling"
        if (params.skip_qc)             summary['QC tools skipped']                  = skip_qc.join(', ')
        if (params.target_bed)          summary['Target BED']                        = params.target_bed
        if (params.tools)               summary['Tools']                             = tools.join(', ')
        if (params.trim_fastq || params.split_fastq) summary['Modify fastqs'] = "trim and/or split"

        if (params.trim_fastq) {
            summary['Fastq trim']         = "Fastq trim selected"
            summary['Trim R1']            = "${params.clip_r1} bp"
            summary['Trim R2']            = "${params.clip_r2} bp"
            summary["Trim 3 R1"]         = "${params.three_prime_clip_r1} bp"
            summary["Trim 3 R2"]         = "${params.three_prime_clip_r2} bp"
            summary['NextSeq Trim']       = "${params.trim_nextseq} bp"
            summary['Saved Trimmed Fastq'] = params.save_trimmed ? 'Yes' : 'No'
        }
        if (params.split_fastq)          summary['Reads in fastq']                   = params.split_fastq

        summary['MarkDuplicates'] = "Options"
        summary['Java options'] = params.markdup_java_options
        summary['GATK Spark']   = params.use_gatk_spark ? 'Yes' : 'No'

        summary['Save BAMs mapped']   = params.save_bam_mapped ? 'Yes' : 'No'
        summary['Skip MarkDuplicates']   = params.skip_markduplicates ? 'Yes' : 'No'

        if ('ascat' in tools) {
            summary['ASCAT'] = "Options"
            if (params.ascat_purity) summary['purity'] = params.ascat_purity
            if (params.ascat_ploidy) summary['ploidy'] = params.ascat_ploidy
        }

        if ('controlfreec' in tools) {
            summary['Control-FREEC'] = "Options"
            if (params.cf_window)    summary['window']             = params.cf_window
            if (params.cf_coeff)     summary['coeff of variation'] = params.cf_coeff
            if (params.cf_ploidy)    summary['ploidy']             = params.cf_ploidy
        }

        if ('haplotypecaller' in tools)             summary['GVCF']       = params.generate_gvcf ? 'Yes' : 'No'
        if ('strelka' in tools && 'manta' in tools) summary['Strelka BP'] = params.no_strelka_bp ? 'No' : 'Yes'
        if (params.pon && ('mutect2' in tools || (params.sentieon && 'tnscope' in tools))) summary['Panel of normals'] = params.pon

        if (params.annotate_tools) summary['Tools to annotate'] = annotate_tools.join(', ')

        if (params.annotation_cache) {
            summary['Annotation cache'] = "Enabled"
            if (params.snpeff_cache) summary['snpEff cache'] = params.snpeff_cache
            if (params.vep_cache)    summary['VEP cache']    = params.vep_cache
        }

        if (params.cadd_cache) {
            summary['CADD cache'] = "Enabled"
            if (params.cadd_indels)  summary['CADD indels']  = params.cadd_indels
            if (params.cadd_wg_snvs) summary['CADD wg snvs'] = params.cadd_wg_snvs
        }

        if (params.genesplicer) summary['genesplicer'] = "Enabled"

        if (params.igenomes_base && !params.igenomes_ignore) summary['AWS iGenomes base'] = params.igenomes_base
        if (params.igenomes_ignore)                          summary['AWS iGenomes']      = "Do not use"
        if (params.genomes_base && !params.igenomes_ignore)  summary['Genomes base']      = params.genomes_base

        summary['Save Reference']    = params.save_reference ? 'Yes' : 'No'

        if (params.ac_loci)                 summary['Loci']                    = params.ac_loci
        if (params.ac_loci_gc)              summary['Loci GC']                 = params.ac_loci_gc
        if (params.bwa)                     summary['BWA indexes']             = params.bwa
        if (params.chr_dir)                 summary['Chromosomes']             = params.chr_dir
        if (params.chr_length)              summary['Chromosomes length']      = params.chr_length
        if (params.dbsnp)                   summary['dbsnp']                   = params.dbsnp
        if (params.dbsnp_index)             summary['dbsnp index']             = params.dbsnp_index
        if (params.dict)                    summary['dict']                    = params.dict
        if (params.fasta)                   summary['fasta reference']         = params.fasta
        if (params.fasta_fai)               summary['fasta index']             = params.fasta_fai
        if (params.germline_resource)       summary['germline resource']       = params.germline_resource
        if (params.germline_resource_index) summary['germline resource index'] = params.germline_resource_index
        if (params.intervals)               summary['intervals']               = params.intervals
        if (params.known_indels)            summary['known indels']            = params.known_indels
        if (params.known_indels_index)      summary['known indels index']      = params.known_indels_index
        if (params.mappability)             summary['Mappability']             = params.mappability
        if (params.snpeff_cache)            summary['snpEff cache']            = params.snpeff_cache
        if (params.snpeff_db)               summary['snpEff DB']               = params.snpeff_db
        if (params.species)                 summary['snpEff species']          = params.species
        if (params.vep_cache)               summary['VEP cache']               = params.vep_cache
        if (params.vep_cache_version)       summary['VEP cache version']       = params.vep_cache_version

        summary['Output dir']       = params.outdir
        summary['Publish dir mode'] = params.publish_dir_mode
        if (params.sequencing_center) summary['Sequenced by'] = params.sequencing_center

        summary['Launch dir']       = workflow.launchDir
        summary['Working dir']      = workflow.workDir
        summary['Script dir']       = workflow.projectDir
        summary['User']             = workflow.userName

        if (params.multiqc_config) summary['MultiQC config'] = params.multiqc_config

        summary['Config Profile'] = workflow.profile

        if (params.config_profile_description) summary['Description']   = params.config_profile_description
        if (params.config_profile_contact)     summary['Contact'] = params.config_profile_contact
        if (params.config_profile_url)         summary['URL']     = params.config_profile_url

        summary['Config Files'] = workflow.configFiles.join(', ')

        if (params.email || params.email_on_fail) {
            summary['E-mail Address']    = params.email
            summary['E-mail on failure'] = params.email_on_fail
            summary['MultiQC maxsize']   = params.max_multiqc_email_size
        }

        if (workflow.profile.contains('awsbatch')) {
            summary['AWS Region']   = params.awsregion
            summary['AWS Queue']    = params.awsqueue
            summary['AWS CLI']      = params.awscli
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
