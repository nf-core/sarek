//
// Run GATK mutect2 in tumor only mode, getepileupsummaries, calculatecontamination and filtermutectcalls
//

params.mutect2_options     = [:]
params.getpileup_options   = [:]
params.calccontam_options  = [:]
params.filtercalls_options = [suffix: '_filtered']

include { GATK4_MUTECT2                as MUTECT2 }                  from '../../../modules/gatk4/mutect2/main'                addParams( options: params.mutect2_options )
include { GATK4_GETPILEUPSUMMARIES     as GETPILEUPSUMMARIES }       from '../../../modules/gatk4/getpileupsummaries/main'     addParams( options: params.getpileup_options )
include { GATK4_CALCULATECONTAMINATION as CALCULATECONTAMINATION }   from '../../../modules/gatk4/calculatecontamination/main' addParams( options: params.calccontam_options )
include { GATK4_FILTERMUTECTCALLS      as FILTERMUTECTCALLS }        from '../../../modules/gatk4/filtermutectcalls/main'      addParams( options: params.filtercalls_options )

workflow GATK_TUMOR_ONLY_SOMATIC_VARIANT_CALLING {
    take:
    input                     // channel: [ val(meta), [ input ], [ input_index ], [] ]
    fasta                     // channel: /path/to/reference/fasta
    fai                       // channel: /path/to/reference/fasta/index
    dict                      // channel: /path/to/reference/fasta/dictionary
    germline_resource         // channel: /path/to/germline/resource
    germline_resource_tbi     // channel: /path/to/germline/index
    panel_of_normals          // channel: /path/to/panel/of/normals
    panel_of_normals_tbi      // channel: /path/to/panel/of/normals/index
    interval_file             // channel: /path/to/interval/file


    main:
    ch_versions = Channel.empty()
    mutect2_input = channel.from(input)

    //
    //Perform variant calling using mutect2 module in tumor single mode.
    //
    MUTECT2 ( mutect2_input , true , false , false , [] , fasta , fai , dict , germline_resource , germline_resource_tbi , panel_of_normals , panel_of_normals_tbi )
    ch_versions = ch_versions.mix(MUTECT2.out.versions)

    //
    //Generate pileup summary table using getepileupsummaries.
    //
    pileup_input = channel.from(input).map {
        meta, input_file, input_index, which_norm ->
        [meta, input_file[0], input_index[0]]
    }
    GETPILEUPSUMMARIES ( pileup_input , germline_resource , germline_resource_tbi , interval_file )
    ch_versions = ch_versions.mix(GETPILEUPSUMMARIES.out.versions)

    //
    //Contamination and segmentation tables created using calculatecontamination on the pileup summary table.
    //
    ch_pileup = GETPILEUPSUMMARIES.out.table.collect()
    //[] is a placeholder for the optional input where the matched normal sample would be passed in for tumor-normal samples, which is not necessary for this workflow.
    ch_pileup.add([])
    CALCULATECONTAMINATION ( ch_pileup, true )
    ch_versions = ch_versions.mix(CALCULATECONTAMINATION.out.versions)

    //
    //Mutect2 calls filtered by filtermutectcalls using the contamination and segmentation tables.
    //
    ch_vcf =           MUTECT2.out.vcf.collect()
    ch_tbi =           MUTECT2.out.tbi.collect()
    ch_stats =         MUTECT2.out.stats.collect()
    //[] is added as a placeholder for the optional input file artifact priors, which is only used for tumor-normal samples and therefor isn't needed in this workflow.
    ch_stats.add([])
    ch_segment =       CALCULATECONTAMINATION.out.segmentation.collect()
    ch_contamination = CALCULATECONTAMINATION.out.contamination.collect()
    //[] is added as a placeholder for entering a contamination estimate value, which is not needed as this workflow uses the contamination table instead.
    ch_contamination.add([])
    ch_filtermutect_in = ch_vcf.combine(ch_tbi, by: 0).combine(ch_stats, by: 0).combine(ch_segment, by: 0).combine(ch_contamination, by: 0)
    FILTERMUTECTCALLS ( ch_filtermutect_in, fasta, fai, dict )
    ch_versions = ch_versions.mix(FILTERMUTECTCALLS.out.versions)

    emit:
    mutect2_vcf            = MUTECT2.out.vcf.collect()                             // channel: [ val(meta), [ vcf ] ]
    mutect2_index          = MUTECT2.out.tbi.collect()                             // channel: [ val(meta), [ tbi ] ]
    mutect2_stats          = MUTECT2.out.stats.collect()                           // channel: [ val(meta), [ stats ] ]

    pileup_table           = GETPILEUPSUMMARIES.out.table.collect()                // channel: [ val(meta), [ table ] ]

    contamination_table    = CALCULATECONTAMINATION.out.contamination.collect()    // channel: [ val(meta), [ contamination ] ]
    segmentation_table     = CALCULATECONTAMINATION.out.segmentation.collect()     // channel: [ val(meta), [ segmentation ] ]

    filtered_vcf           = FILTERMUTECTCALLS.out.vcf.collect()                   // channel: [ val(meta), [ vcf ] ]
    filtered_index         = FILTERMUTECTCALLS.out.tbi.collect()                   // channel: [ val(meta), [ tbi ] ]
    filtered_stats         = FILTERMUTECTCALLS.out.stats.collect()                 // channel: [ val(meta), [ stats ] ]

    versions               = ch_versions                                           // channel: [ versions.yml ]
}
