include { LOFREQ_SOMATIC  as  LOFREQSOMATIC } from '../../../modules/nf-core/lofreq/somatic/main.nf'
include { LOFREQ_SOMATIC  as  LOFREQSOMATIC_FILTERED } from '../../../modules/nf-core/lofreq/somatic/main.nf'
include { GATK4_MERGEVCFS as MERGE_LOFREQ_INDELS } from '../../../modules/nf-core/gatk4/mergevcfs/main'
include { GATK4_MERGEVCFS as MERGE_LOFREQ_SNVS   } from '../../../modules/nf-core/gatk4/mergevcfs/main'
include { BCFTOOLS_FILTER as FILTER } from '../../../modules/nf-core/bcftools/filter/main'

workflow BAM_VARIANT_CALLING_SOMATIC_LOFREQ {
    take:
    input         // channel: [mandatory] [ meta, normal_cram, normal_crai, tumor_cram, tumor_crai ]
    dict          // channel: [optional]  [ meta, dict ]
    fasta         // channel: [mandatory] [ fasta ]
    fai           // channel: [mandatory] [ fasta_fai ]
    dbsnp         // channel: /path/to/dbsnp/resource
    dbsnp_tbi     // channel: /path/to/dbsnp/index
    intervals     // channel: [mandatory] [ intervals ]

    main:
    versions = Channel.empty()
    


    input_intervals = input.combine(intervals)
        // Sort channel elements for LOFREQ_SOMATIC module
        .map {meta, normal_cram, normal_crai, tumor_cram, tumor_crai, intervals, num_intervals -> [meta + [num_intervals:num_intervals], tumor_cram, tumor_crai, normal_cram, normal_crai, intervals]}

        

    if (params.dbsnp && params.dbsnp_lofreq){

        // ===================================================================================================
        // Organizar los inputs de FILTER
        // ===================================================================================================
        input_filter = dbsnp ? dbsnp.map{ dbsnp -> [ [ id:dbsnp.baseName ], dbsnp ] } : [[id: 'null'], []]
        input_filter.subscribe onNext: { println it }, onComplete: { println 'Done' }
        FILTER(input_filter)

        // ===================================================================================================
        // Administrar outputs de FILTER
        // =================================================================================================== 
        dbsnp_new = FILTER.out.vcf
        dbsnp_new_tbi = FILTER.out.tbi
        

        // ===================================================================================================
        // Printear los parámetros 
        // ===================================================================================================
        dbsnp_new.subscribe onNext: { println it }, onComplete: { println 'Done' } 
        dbsnp_new_tbi.subscribe onNext: { println it }, onComplete: { println 'Done' }
        input_intervals.subscribe onNext: { println it }, onComplete: { println 'Done' }
        fasta.subscribe onNext: { println it }, onComplete: { println 'Done' }
        fai.subscribe onNext: { println it }, onComplete: { println 'Done' }

        // ===================================================================================================
        // Ejecutar lofreq somátic
        // ===================================================================================================
        LOFREQSOMATIC_FILTERED(input_intervals, fasta, fai, dbsnp_new, dbsnp_new_tbi)
        final_indels = LOFREQSOMATIC_FILTERED.out.vcf_indels
        final_snvs = LOFREQSOMATIC_FILTERED.out.vcf_snvs
        final_version = LOFREQSOMATIC_FILTERED.out.versions

    }else{
        LOFREQSOMATIC(input_intervals, fasta, fai, '', '')
        final_indels = LOFREQSOMATIC.out.vcf_indels
        final_snvs = LOFREQSOMATIC.out.vcf_snvs
        final_version = LOFREQSOMATIC.out.versions
    }

    // Figuring out if there is one or more vcf(s) from the same sample
    vcf_indels = final_indels.branch{
        // Use meta.num_intervals to asses number of intervals
        intervals:    it[0].num_intervals > 1
        no_intervals: it[0].num_intervals <= 1
    }

    // Figuring out if there is one or more vcf(s) from the same sample
    vcf_snvs = final_snvs.branch{
        // Use meta.num_intervals to asses number of intervals
        intervals:    it[0].num_intervals > 1
        no_intervals: it[0].num_intervals <= 1
    }

    // Only when using intervals
    vcf_indels_to_merge = vcf_indels.intervals.map{ meta, vcf -> [ groupKey(meta, meta.num_intervals), vcf ]}.groupTuple()
    vcf_snvs_to_merge = vcf_snvs.intervals.map{ meta, vcf -> [ groupKey(meta, meta.num_intervals), vcf ]}.groupTuple()

    MERGE_LOFREQ_INDELS(vcf_indels_to_merge, dict)
    MERGE_LOFREQ_SNVS(vcf_snvs_to_merge, dict)

    // Mix intervals and no_intervals channels together
    vcf = Channel.empty().mix(MERGE_LOFREQ_INDELS.out.vcf, MERGE_LOFREQ_SNVS.out.vcf, vcf_indels.no_intervals, vcf_snvs.no_intervals)
        // add variantcaller to meta map and remove no longer necessary field: num_intervals
        .map{ meta, vcf -> [ meta - meta.subMap('num_intervals') + [ variantcaller:'lofreq' ], vcf ] }

    versions = versions.mix(final_version)

    emit:
    vcf
    versions
}
