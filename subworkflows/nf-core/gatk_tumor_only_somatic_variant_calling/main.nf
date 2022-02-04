//
// Run GATK mutect2 in tumor only mode, getepileupsummaries, calculatecontamination and filtermutectcalls
//

include { BGZIP as BGZIP_MUTECT2                      } from '../../../modules/local/bgzip'
include { CONCAT_VCF as CONCAT_VCF_MUTECT2            } from '../../../modules/local/concat_vcf/main'

include { GATK4_MUTECT2                as MUTECT2 }                  from '../../../modules/nf-core/modules/gatk4/mutect2/main'
include { GATK4_MERGEMUTECTSTATS       as MERGEMUTECTSTATS }         from '../../../modules/local/gatk4/mergemutectstats'
include { GATK4_GETPILEUPSUMMARIES     as GETPILEUPSUMMARIES }       from '../../../modules/nf-core/modules/gatk4/getpileupsummaries/main'
include { GATK4_GATHERPILEUPSUMMARIES  as GATHERPILEUPSUMMARIES }    from '../../../modules/local/gatk4/gatherpileupsummaries'
include { GATK4_CALCULATECONTAMINATION as CALCULATECONTAMINATION }   from '../../../modules/nf-core/modules/gatk4/calculatecontamination/main'
include { GATK4_FILTERMUTECTCALLS      as FILTERMUTECTCALLS }        from '../../../modules/nf-core/modules/gatk4/filtermutectcalls/main'

workflow GATK_TUMOR_ONLY_SOMATIC_VARIANT_CALLING {
    take:
    input                     // channel: [ val(meta), [ input ], [ input_index ], [intervals], [] ]
    fasta                     // channel: /path/to/reference/fasta
    fai                       // channel: /path/to/reference/fasta/index
    dict                      // channel: /path/to/reference/fasta/dictionary
    germline_resource         // channel: /path/to/germline/resource
    germline_resource_tbi     // channel: /path/to/germline/index
    panel_of_normals          // channel: /path/to/panel/of/normals
    panel_of_normals_tbi      // channel: /path/to/panel/of/normals/index
    num_intervals
    no_intervals
    intervals_bed_combine_gz


    main:
    ch_versions = Channel.empty()

    //
    //Perform variant calling using mutect2 module in tumor single mode.
    //
    MUTECT2 ( input , true , false , false , [] , fasta , fai , dict , germline_resource , germline_resource_tbi , panel_of_normals , panel_of_normals_tbi )
    ch_versions = ch_versions.mix(MUTECT2.out.versions)

    //
    //Generate pileup summary table using getepileupsummaries.
    //
    pileup_input = input.map {
        meta, input_file, input_index, intervals, which_norm ->
        [meta, input_file, input_index, intervals]
    }
    GETPILEUPSUMMARIES ( pileup_input , fasta, fai, dict, germline_resource , germline_resource_tbi )
    ch_versions = ch_versions.mix(GETPILEUPSUMMARIES.out.versions)

    if(no_intervals){
        mutect2_vcf_gz_tbi = MUTECT2.out.vcf.join(MUTECT2.out.tbi)
        mutect2_stats = MUTECT2.out.stats
        pileup_table = GETPILEUPSUMMARIES.out.table
    }else{

        //Merge Mutect2 VCF
        BGZIP_MUTECT2(MUTECT2.out.vcf)
        BGZIP_MUTECT2.out.vcf.map{ meta, vcf ->
            new_meta = meta.clone()
            new_meta.id = new_meta.sample
            [new_meta, vcf]
        }.set{bgzip_mutect2}

        mutect2_vcf_to_concat = bgzip_mutect2.groupTuple(size: num_intervals)

        CONCAT_VCF_MUTECT2(mutect2_vcf_to_concat, fai, intervals_bed_combine_gz)
        mutect2_vcf_gz_tbi = CONCAT_VCF_MUTECT2.out.vcf

        ch_versions = ch_versions.mix(BGZIP_MUTECT2.out.versions)
        ch_versions = ch_versions.mix(CONCAT_VCF_MUTECT2.out.versions)

        //Merge Muteect2 Stats
        MUTECT2.out.stats.map{ meta, stats ->
            new_meta = meta.clone()
            new_meta.id = new_meta.sample
            [new_meta, stats]
        }.groupTuple(size: num_intervals).set{mutect2_stats_to_merge}

        MERGEMUTECTSTATS(mutect2_stats_to_merge)
        mutect2_stats = MERGEMUTECTSTATS.out.stats
        ch_versions = ch_versions.mix(MERGEMUTECTSTATS.out.versions)

        //Merge Pileup Summaries
        pileup_tables_to_gather = GETPILEUPSUMMARIES.out.table.map{ meta, table ->
            new_meta = meta.clone()
            new_meta.id = new_meta.sample
            [new_meta, table]
        }.groupTuple(size: num_intervals)

        pileup_tables_to_gather.view()

        GATHERPILEUPSUMMARIES(pileup_tables_to_gather, dict)
        pileup_table = GATHERPILEUPSUMMARIES.out.table

    }

    //
    //Contamination and segmentation tables created using calculatecontamination on the pileup summary table.
    //
    pileup_table.map{meta, table -> [meta, table, []]}.set{table_contamination}
    CALCULATECONTAMINATION ( table_contamination, true )
    ch_versions = ch_versions.mix(CALCULATECONTAMINATION.out.versions)


    //
    //Mutect2 calls filtered by filtermutectcalls using the contamination and segmentation tables.
    //
    ch_filtermutect = mutect2_vcf_gz_tbi.join(mutect2_stats)
                                            .join(CALCULATECONTAMINATION.out.segmentation)
                                            .join(CALCULATECONTAMINATION.out.contamination)
    ch_filtermutect.view()
    ch_filtermutect_in = ch_filtermutect.map{ meta, vcf, tbi, stats, seg, cont -> [meta, vcf, tbi, stats, [], seg, cont, []] }
    FILTERMUTECTCALLS ( ch_filtermutect_in, fasta, fai, dict )
    ch_versions = ch_versions.mix(FILTERMUTECTCALLS.out.versions)

    emit:
    mutect2_vcf_gz_tbi            = mutect2_vcf_gz_tbi               // channel: [ val(meta), [ vcf ] ]
    mutect2_stats                 = MUTECT2.out.stats                // channel: [ val(meta), [ stats ] ]

    pileup_table                  = pileup_table                     // channel: [ val(meta), [ table ] ]

    contamination_table    = CALCULATECONTAMINATION.out.contamination    // channel: [ val(meta), [ contamination ] ]
    segmentation_table     = CALCULATECONTAMINATION.out.segmentation     // channel: [ val(meta), [ segmentation ] ]

    filtered_vcf           = FILTERMUTECTCALLS.out.vcf                   // channel: [ val(meta), [ vcf ] ]
    filtered_index         = FILTERMUTECTCALLS.out.tbi                   // channel: [ val(meta), [ tbi ] ]
    filtered_stats         = FILTERMUTECTCALLS.out.stats                 // channel: [ val(meta), [ stats ] ]

    versions               = ch_versions                                           // channel: [ versions.yml ]
}
