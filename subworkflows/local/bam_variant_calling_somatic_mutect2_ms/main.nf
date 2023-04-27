//
// Run GATK mutect2 in tumor normal mode, getepileupsummaries, calculatecontamination, learnreadorientationmodel and filtermutectcalls
//

workflow BAM_VARIANT_CALLING_SOMATIC_MUTECT2_MS {
    take:
    input                     // channel: [ val(meta), [ input ], [ input_index ] ]
    fasta                     // channel: /path/to/reference/fasta
    fai                       // channel: /path/to/reference/fasta/index
    dict                      // channel: /path/to/reference/fasta/dictionary
    germline_resource         // channel: /path/to/germline/resource
    germline_resource_tbi     // channel: /path/to/germline/index
    panel_of_normals          // channel: /path/to/panel/of/normals
    panel_of_normals_tbi      // channel: /path/to/panel/of/normals/index
    intervals                 // channel: [mandatory] intervals/target regions
    intervals_bed_gz_tbi      // channel: [mandatory] intervals/target regions index zipped and indexed
    intervals_bed_combined    // channel: [mandatory] intervals/target regions in one file unzipped

    main:
    ch_versions = Channel.empty()

    // Multi sample somatic variant calling channels
    // Extract all cram files for each patient
    input.branch{
            normal: it[0].status == 0
            tumor:  it[0].status == 1
        }.set{ch_cram_variant_calling_branched}
    ch_cram_variant_calling_to_cross = input.map{ meta, cram, crai -> [meta.patient, cram, crai] }.groupTuple()

    // Collect tumor and normal sample names separately
    ch_cram_variant_calling_normal_names = ch_cram_variant_calling_branched.normal.map{ meta, cram, crai-> [meta.patient, meta.id] }.groupTuple()
    ch_cram_variant_calling_tumor_names = ch_cram_variant_calling_branched.tumor.map{ meta, cram, crai-> [meta.patient, meta.id] }.groupTuple()
    
    // Join sample name channels with cram channels by patient
    ch_cram_variant_calling_with_normal_names = ch_cram_variant_calling_to_cross.join(ch_cram_variant_calling_normal_names, remainder: true)
    ch_cram_variant_calling_for_multi_sample_mutect2 = ch_cram_variant_calling_with_normal_names.join(ch_cram_variant_calling_tumor_names, remainder: true)
    //ch_cram_variant_calling_for_multi_sample_mutect2.view()

    // Separate patients with tumor sample(s)
    ch_cram_variant_calling_for_multi_sample_mutect2.branch{
        somatic: it[4] != null
        other: it[4] == null
    }.set{ch_cram_variant_calling_for_multi_sample_mutect2_branched}
    ch_cram_variant_calling_for_multi_sample_mutect2_branched.somatic.view()

    emit:
    versions               = ch_versions                                    // channel: [ versions.yml ]
}
