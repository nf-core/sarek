//
// SENTIEON HAPLOTYPER germline variant calling
//
// For all modules here:
// A when clause condition is defined in the conf/modules.config to determine if the module should be run

include { GATK4_MERGEVCFS            as MERGE_SENTIEON_HAPLOTYPER_GVCFS } from '../../../modules/nf-core/gatk4/mergevcfs/main'
include { GATK4_MERGEVCFS            as MERGE_SENTIEON_HAPLOTYPER_VCFS  } from '../../../modules/nf-core/gatk4/mergevcfs/main'
include { SENTIEON_HAPLOTYPER                                           } from '../../../modules/nf-core/sentieon/haplotyper/main'

workflow BAM_VARIANT_CALLING_SENTIEON_HAPLOTYPER {
    take:
    cram                           // channel: [mandatory] [ meta, cram, crai, interval.bed ]
    fasta                          // channel: [mandatory]
    fasta_fai                      // channel: [mandatory]
    dict                           // channel: [mandatory]
    dbsnp                          // channel: [optional]
    dbsnp_tbi                      // channel: [optional]
    dbsnp_vqsr                     // channel: [optional]
    intervals                      // channel: [mandatory] [ intervals, num_intervals ] or [ [], 0 ] if no intervals
    joint_germline                 // boolean: [mandatory] [default: false] joint calling of germline variants
    sentieon_haplotyper_emit_mode

    main:
    versions = Channel.empty()

    gvcf               = Channel.empty()
    vcf                = Channel.empty()
    genotype_intervals = Channel.empty()

    // Combine cram and intervals for spread and gather strategy
    cram_intervals_for_sentieon = cram.combine(intervals)
        // Move num_intervals to meta map
        .map{ meta, cram, crai, intervals, num_intervals -> [
            meta + [
                num_intervals:num_intervals,
                intervals_name:intervals.baseName,
                variantcaller:'sentieon_haplotyper'],
            cram,
            crai,
            intervals
            ]
        }


    emit_mode_items = sentieon_haplotyper_emit_mode.split(',').each{ it -> it.toLowerCase().trim() }
    lst = emit_mode_items - 'gvcf'
    emit_vcf = lst.size() > 0 ? lst[0] : ''

    SENTIEON_HAPLOTYPER(
        cram_intervals_for_sentieon.map{ meta, cram, crai, intervals, num_intervals -> [ meta, cram, crai, intervals, [] ]},
        fasta,
        fasta_fai,
        dbsnp,
        dbsnp_tbi,
        emit_vcf,
        emit_mode_items.any{ it.equals('gvcf') })

    if (joint_germline) {
        genotype_intervals = SENTIEON_HAPLOTYPER.out.gvcf
            .join(SENTIEON_HAPLOTYPER.out.gvcf_tbi, failOnMismatch: true)
            .join(cram_intervals_for_sentieon, failOnMismatch: true)
            .map{ meta, gvcf, tbi, cram, crai, intervals -> [ meta, gvcf, tbi, intervals ] }
    }

    // Figure out if using intervals or no_intervals
    haplotyper_vcf_branch = SENTIEON_HAPLOTYPER.out.vcf.map{
            meta, vcf -> [ meta - meta.subMap('interval_name'), vcf]
        }
        .branch{
            intervals:    it[0].num_intervals > 1
            no_intervals: it[0].num_intervals <= 1
        }

    haplotyper_vcf_tbi_branch = SENTIEON_HAPLOTYPER.out.vcf_tbi.map{
            meta, vcf_tbi -> [ meta - meta.subMap('interval_name'), vcf_tbi]
        }
        .branch{
            intervals:    it[0].num_intervals > 1
            no_intervals: it[0].num_intervals <= 1
        }

    haplotyper_gvcf_branch = SENTIEON_HAPLOTYPER.out.gvcf.map{
            meta, gvcf -> [ meta - meta.subMap('interval_name'), gvcf]
        }
        .branch{
            intervals:    it[0].num_intervals > 1
            no_intervals: it[0].num_intervals <= 1
        }

    haplotyper_gvcf_tbi_branch = SENTIEON_HAPLOTYPER.out.gvcf_tbi.map{
            meta, gvcf_tbi -> [ meta - meta.subMap('interval_name'), gvcf_tbi]
        }
        .branch{
            intervals:    it[0].num_intervals > 1
            no_intervals: it[0].num_intervals <= 1
        }

    vcfs_for_merging = haplotyper_vcf_branch.intervals.map{
        meta, vcf -> [ groupKey(meta, meta.num_intervals), vcf ]}

    vcfs_for_merging = vcfs_for_merging.map{
        meta, vcf -> [
            meta - meta.subMap('intervals_name'),
            vcf]}.groupTuple()

    // VCFs
    // Only when using intervals
    MERGE_SENTIEON_HAPLOTYPER_VCFS(vcfs_for_merging, dict)

    haplotyper_vcf = Channel.empty().mix(
        MERGE_SENTIEON_HAPLOTYPER_VCFS.out.vcf,
        haplotyper_vcf_branch.no_intervals)

    haplotyper_tbi = Channel.empty().mix(
        MERGE_SENTIEON_HAPLOTYPER_VCFS.out.tbi,
        haplotyper_vcf_tbi_branch.no_intervals)

    // Remove no longer necessary field: num_intervals
    vcf = haplotyper_vcf.map{ meta, vcf -> [ meta - meta.subMap('num_intervals'), vcf ] }
    vcf_tbi = haplotyper_tbi.map{ meta, tbi -> [ meta - meta.subMap('num_intervals'), tbi ] }

    // GVFs
    // Only when using intervals
    gvcfs_for_merging = haplotyper_gvcf_branch.intervals.map{
        meta, vcf -> [groupKey(meta, meta.num_intervals), vcf]}

    gvcfs_for_merging = gvcfs_for_merging.map{
        meta, vcf -> [ meta - meta.subMap('intervals_name'), vcf ]
    }.groupTuple()

    MERGE_SENTIEON_HAPLOTYPER_GVCFS(gvcfs_for_merging, dict)

    gvcf = Channel.empty().mix(
        MERGE_SENTIEON_HAPLOTYPER_GVCFS.out.vcf,
        haplotyper_gvcf_branch.no_intervals)

    gvcf_tbi = Channel.empty().mix(
        MERGE_SENTIEON_HAPLOTYPER_GVCFS.out.tbi,
        haplotyper_gvcf_tbi_branch.no_intervals)

    versions = versions.mix(SENTIEON_HAPLOTYPER.out.versions)
    versions = versions.mix(MERGE_SENTIEON_HAPLOTYPER_VCFS.out.versions)
    versions = versions.mix(MERGE_SENTIEON_HAPLOTYPER_GVCFS.out.versions)

    emit:
    versions
    vcf
    vcf_tbi
    gvcf
    gvcf_tbi
    genotype_intervals // For joint genotyping

}
