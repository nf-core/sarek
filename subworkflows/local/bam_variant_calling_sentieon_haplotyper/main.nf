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
        .map{ meta, cram_, crai, intervals_, num_intervals -> [
            meta + [
                num_intervals:num_intervals,
                intervals_name:intervals_.baseName,
                variantcaller:'sentieon_haplotyper'],
            cram_,
            crai,
            intervals_
            ]
        }

    emit_mode_items = sentieon_haplotyper_emit_mode.split(',').each{ it -> it.toLowerCase().trim() }
    lst = emit_mode_items - 'gvcf'
    emit_vcf = lst.size() > 0 ? lst[0] : ''

    SENTIEON_HAPLOTYPER(
        cram_intervals_for_sentieon.map{ meta, cram_, crai, intervals_ -> [ meta, cram_, crai, intervals_, [] ]},
        fasta,
        fasta_fai,
        dbsnp.map{file -> [[id:'dbsnp'], file]},
        dbsnp_tbi.map{file -> [[id:'dbsnp'], file]},
        emit_vcf,
        emit_mode_items.any{ it.equals('gvcf') })

    if (joint_germline) {
        genotype_intervals = SENTIEON_HAPLOTYPER.out.gvcf
            .join(SENTIEON_HAPLOTYPER.out.gvcf_tbi, failOnMismatch: true)
            .join(cram_intervals_for_sentieon, failOnMismatch: true)
            .map{ meta, gvcf_, tbi, _cram, _crai, intervals_ -> [ meta, gvcf_, tbi, intervals_ ] }
    }

    // Figure out if using intervals or no_intervals.
    // Strip `intervals_name` here so the branch output meta is the same one we'll
    // group on later — fixes the prior `interval_name` typo which was a no-op.
    haplotyper_vcf_branch = SENTIEON_HAPLOTYPER.out.vcf.map{
            meta, vcf_ -> [ meta - meta.subMap('intervals_name'), vcf_]
        }
        .branch{
            intervals:    it[0].num_intervals > 1
            no_intervals: it[0].num_intervals <= 1
        }

    haplotyper_vcf_tbi_branch = SENTIEON_HAPLOTYPER.out.vcf_tbi.map{
            meta, vcf_tbi -> [ meta - meta.subMap('intervals_name'), vcf_tbi]
        }
        .branch{
            intervals:    it[0].num_intervals > 1
            no_intervals: it[0].num_intervals <= 1
        }

    haplotyper_gvcf_branch = SENTIEON_HAPLOTYPER.out.gvcf.map{
            meta, gvcf_ -> [ meta - meta.subMap('intervals_name'), gvcf_]
        }
        .branch{
            intervals:    it[0].num_intervals > 1
            no_intervals: it[0].num_intervals <= 1
        }

    haplotyper_gvcf_tbi_branch = SENTIEON_HAPLOTYPER.out.gvcf_tbi.map{
            meta, gvcf_tbi -> [ meta - meta.subMap('intervals_name'), gvcf_tbi]
        }
        .branch{
            intervals:    it[0].num_intervals > 1
            no_intervals: it[0].num_intervals <= 1
        }

    // Per-sample merge. Wrap the (already-`intervals_name`-stripped) meta in
    // groupKey *inside* the map so the GroupKey survives into groupTuple —
    // `groupKey` followed by a separate `meta - subMap(...)` step would unwrap
    // it (Map.minus(Map) returns a plain LinkedHashMap), losing the size hint
    // and forcing groupTuple to wait until the upstream channel closes. With
    // it preserved, per-sample merges emit progressively.
    vcfs_for_merging = haplotyper_vcf_branch.intervals.map{
        meta, vcf_ -> [ groupKey(meta, meta.num_intervals), vcf_ ]
    }.groupTuple()

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
    vcf = haplotyper_vcf.map{ meta, vcf_ -> [ meta - meta.subMap('num_intervals'), vcf_ ] }
    vcf_tbi = haplotyper_tbi.map{ meta, tbi -> [ meta - meta.subMap('num_intervals'), tbi ] }

    // GVCFs — see vcfs_for_merging comment above for why the groupKey is built
    // inside the same map (preserves size hint into groupTuple).
    gvcfs_for_merging = haplotyper_gvcf_branch.intervals.map{
        meta, vcf_ -> [ groupKey(meta, meta.num_intervals), vcf_ ]
    }.groupTuple()

    MERGE_SENTIEON_HAPLOTYPER_GVCFS(gvcfs_for_merging, dict)

    gvcf = Channel.empty().mix(
        MERGE_SENTIEON_HAPLOTYPER_GVCFS.out.vcf,
        haplotyper_gvcf_branch.no_intervals)

    gvcf_tbi = Channel.empty().mix(
        MERGE_SENTIEON_HAPLOTYPER_GVCFS.out.tbi,
        haplotyper_gvcf_tbi_branch.no_intervals)


    emit:
    versions
    vcf
    vcf_tbi
    gvcf
    gvcf_tbi
    genotype_intervals // For joint genotyping

}
