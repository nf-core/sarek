//
// DOWNLOAD CACHE SNPEFF VEP
//

include { MULTIQC                 } from '../../../modules/nf-core/multiqc/main'
include { paramsSummaryMap        } from 'plugin/nf-validation'
include { methodsDescriptionText  } from '../utils_nfcore_sarek_pipeline'
include { paramsSummaryMultiqc    } from '../../nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML  } from '../../nf-core/utils_nfcore_pipeline'

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
