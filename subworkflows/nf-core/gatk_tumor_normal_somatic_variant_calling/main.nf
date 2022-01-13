//
// Run GATK mutect2 in tumor normal mode, getepileupsummaries, calculatecontamination, learnreadorientationmodel and filtermutectcalls
//

params.mutect2_options          = [:]
params.learnorientation_options = [:]
params.getpileup_tumor_options  = [suffix: '_tumor']
params.getpileup_normal_options = [suffix: '_normal']
params.calccontam_options       = [:]
params.filtercalls_options      = [suffix: '_filtered']

include { GATK4_MUTECT2                   as MUTECT2 }                   from '../../../modules/gatk4/mutect2/main'                   addParams( options: params.mutect2_options )
include { GATK4_LEARNREADORIENTATIONMODEL as LEARNREADORIENTATIONMODEL } from '../../../modules/gatk4/learnreadorientationmodel/main' addParams( options: params.learnorientation_options )
include { GATK4_GETPILEUPSUMMARIES        as GETPILEUPSUMMARIES_TUMOR }  from '../../../modules/gatk4/getpileupsummaries/main'        addParams( options: params.getpileup_tumor_options )
include { GATK4_GETPILEUPSUMMARIES        as GETPILEUPSUMMARIES_NORMAL}  from '../../../modules/gatk4/getpileupsummaries/main'        addParams( options: params.getpileup_normal_options )
include { GATK4_CALCULATECONTAMINATION    as CALCULATECONTAMINATION }    from '../../../modules/gatk4/calculatecontamination/main'    addParams( options: params.calccontam_options )
include { GATK4_FILTERMUTECTCALLS         as FILTERMUTECTCALLS }         from '../../../modules/gatk4/filtermutectcalls/main'         addParams( options: params.filtercalls_options )

workflow GATK_TUMOR_NORMAL_SOMATIC_VARIANT_CALLING {
    take:
    input                     // channel: [ val(meta), [ input ], [ input_index ], [which_norm] ]
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

    //
    //Perform variant calling using mutect2 module in tumor single mode.
    //
    mutect2_input = channel.from(input)
    MUTECT2 ( mutect2_input, false, false, false, [], fasta, fai, dict, germline_resource, germline_resource_tbi, panel_of_normals, panel_of_normals_tbi )
    ch_versions = ch_versions.mix(MUTECT2.out.versions)

    //
    //Generate artifactpriors using learnreadorientationmodel on the f1r2 output of mutect2.
    //
    ch_learnread_in = MUTECT2.out.f1r2.collect()
    LEARNREADORIENTATIONMODEL (ch_learnread_in)
    ch_versions = ch_versions.mix(LEARNREADORIENTATIONMODEL.out.versions)

    //
    //Generate pileup summary tables using getepileupsummaries. tumor sample should always be passed in as the first input and input list entries of ch_mutect2_in,
    //to ensure correct file order for calculatecontamination.
    //
    pileup_tumor_input = channel.from(input).map {
        meta, input_file, input_index, which_norm ->
        [meta, input_file[0], input_index[0]]
    }

    pileup_normal_input = channel.from(input).map {
        meta, input_file, input_index, which_norm ->
        [meta, input_file[1], input_index[1]]
    }
    GETPILEUPSUMMARIES_TUMOR ( pileup_tumor_input, germline_resource, germline_resource_tbi, interval_file )
    GETPILEUPSUMMARIES_NORMAL ( pileup_normal_input, germline_resource, germline_resource_tbi, interval_file )
    ch_versions = ch_versions.mix(GETPILEUPSUMMARIES_NORMAL.out.versions)

    //
    //Contamination and segmentation tables created using calculatecontamination on the pileup summary table.
    //
    ch_pileup_tumor     = GETPILEUPSUMMARIES_TUMOR.out.table.collect()
    ch_pileup_normal    = GETPILEUPSUMMARIES_NORMAL.out.table.collect()
    ch_calccon_in       = ch_pileup_tumor.combine(ch_pileup_normal, by: 0)
    CALCULATECONTAMINATION ( ch_calccon_in, true )
    ch_versions   = ch_versions.mix(CALCULATECONTAMINATION.out.versions)

    //
    //Mutect2 calls filtered by filtermutectcalls using the artifactpriors, contamination and segmentation tables.
    //
    ch_vcf                   = MUTECT2.out.vcf.collect()
    ch_tbi                   = MUTECT2.out.tbi.collect()
    ch_stats                 = MUTECT2.out.stats.collect()
    ch_orientation           = LEARNREADORIENTATIONMODEL.out.artifactprior.collect()
    ch_segment               = CALCULATECONTAMINATION.out.segmentation.collect()
    ch_contamination         = CALCULATECONTAMINATION.out.contamination.collect()
    //[] is used as a placeholder for optional input to specify the contamination estimate as a value, since the contamination table is used, this is not needed.
    ch_contamination.add([])
    ch_filtermutect_in       = ch_vcf.combine(ch_tbi, by: 0).combine(ch_stats, by: 0).combine(ch_orientation, by: 0).combine(ch_segment, by: 0).combine(ch_contamination, by: 0)
    FILTERMUTECTCALLS ( ch_filtermutect_in, fasta, fai, dict )
    ch_versions              = ch_versions.mix(FILTERMUTECTCALLS.out.versions)

    emit:
    mutect2_vcf            = MUTECT2.out.vcf.collect()                             // channel: [ val(meta), [ vcf ] ]
    mutect2_tbi            = MUTECT2.out.tbi.collect()                             // channel: [ val(meta), [ tbi ] ]
    mutect2_stats          = MUTECT2.out.stats.collect()                           // channel: [ val(meta), [ stats ] ]
    mutect2_f1r2           = MUTECT2.out.f1r2.collect()                            // channel: [ val(meta), [ f1r2 ] ]

    artifact_priors        = LEARNREADORIENTATIONMODEL.out.artifactprior.collect() // channel: [ val(meta), [ artifactprior ] ]

    pileup_table_tumor     = GETPILEUPSUMMARIES_TUMOR.out.table.collect()          // channel: [ val(meta), [ table_tumor ] ]
    pileup_table_normal    = GETPILEUPSUMMARIES_NORMAL.out.table.collect()         // channel: [ val(meta), [ table_normal ] ]

    contamination_table    = CALCULATECONTAMINATION.out.contamination.collect()    // channel: [ val(meta), [ contamination ] ]
    segmentation_table     = CALCULATECONTAMINATION.out.segmentation.collect()     // channel: [ val(meta), [ segmentation ] ]

    filtered_vcf           = FILTERMUTECTCALLS.out.vcf.collect()                   // channel: [ val(meta), [ vcf ] ]
    filtered_tbi           = FILTERMUTECTCALLS.out.tbi.collect()                   // channel: [ val(meta), [ tbi ] ]
    filtered_stats         = FILTERMUTECTCALLS.out.stats.collect()                 // channel: [ val(meta), [ stats ] ]

    versions               = ch_versions                                           // channel: [ versions.yml ]
}
