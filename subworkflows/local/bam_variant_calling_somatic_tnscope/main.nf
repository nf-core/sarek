//
//
// SENTIEON TNSCOPE: tumor-normal mode variantcalling
//

include { SENTIEON_TNSCOPE } from '../../../modules/nf-core/sentieon/tnscope/main'

workflow BAM_VARIANT_CALLING_SOMATIC_TNSCOPE {
    take:
    input                     // channel: [ meta, [ input ], [ input_index ] ]
    fasta                     // channel: /path/to/reference/fasta
    fai                       // channel: /path/to/reference/fasta/index
    germline_resource         // channel: /path/to/germline/resource
    germline_resource_tbi     // channel: /path/to/germline/index
    panel_of_normals          // channel: /path/to/panel/of/normals
    panel_of_normals_tbi      // channel: /path/to/panel/of/normals/index

    main:
    versions = Channel.empty()

    //If no germline resource is provided, then create an empty channel to avoid GetPileupsummaries from being run
    //germline_resource_pileup     = (germline_resource && germline_resource_tbi) ? germline_resource : Channel.empty()
    //germline_resource_pileup_tbi = germline_resource_tbi ?: Channel.empty()

    // Separate normal cram files
    // Extract tumor cram files
    ch_cram = input.multiMap{ meta, cram, crai ->
            normal: [ meta - meta.subMap('tumor_id') , cram[0], crai[0] ]
            tumor:  [ meta - meta.subMap('tumor_id') , cram[1], crai[1] ]
        }

    // Remove duplicates from normal channel and merge normal and tumor crams by patient
    ch_tn_cram =  ch_cram.normal.unique().mix(ch_cram.tumor).groupTuple()

    SENTIEON_TNSCOPE(
        ch_tn_cram,
        fasta.map{fasta_ -> [[id: "fasta"], fasta_]},
        fai.map{fai_ -> [[id: "fasta"], fai_]},
        germline_resource,
        germline_resource_tbi,
        panel_of_normals,
        panel_of_normals_tbi,
        [[],[]], // cosmic
        [[],[]] // cosmic_tbi
    )

    versions = versions.mix(SENTIEON_TNSCOPE.out.versions)

    emit:
    vcf   = SENTIEON_TNSCOPE.out.vcf    // channel: [ meta, vcf ]
    index = SENTIEON_TNSCOPE.out.index  // channel: [ meta, index ]

    versions // channel: [ versions.yml ]
}
