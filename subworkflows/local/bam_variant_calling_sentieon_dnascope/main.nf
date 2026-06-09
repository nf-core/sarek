//
// SENTIEON HAPLOTYPER germline variant calling
//
// For all modules here:
// A when clause condition is defined in the conf/modules.config to determine if the module should be run

include { GATK4_MERGEVCFS            as MERGE_SENTIEON_DNASCOPE_GVCFS } from '../../../modules/nf-core/gatk4/mergevcfs/main'
include { GATK4_MERGEVCFS            as MERGE_SENTIEON_DNASCOPE_VCFS  } from '../../../modules/nf-core/gatk4/mergevcfs/main'
include { SENTIEON_DNASCOPE                                           } from '../../../modules/nf-core/sentieon/dnascope/main'

workflow BAM_VARIANT_CALLING_SENTIEON_DNASCOPE {
    take:
    cram                              // channel: [mandatory] [ meta, cram, crai, interval.bed ]
    fasta                             // channel: [mandatory]
    fasta_fai                         // channel: [mandatory]
    dict                              // channel: [mandatory]
    dbsnp                             // channel: [optional]
    dbsnp_tbi                         // channel: [optional]
    dbsnp_vqsr                        // channel: [optional]
    intervals                         // channel: [mandatory] [ intervals, num_intervals ] or [ [], 0 ] if no intervals
    joint_germline                    // boolean: [mandatory] [default: false] joint calling of germline variants
    sentieon_dnascope_emit_mode       // string
    sentieon_dnascope_pcr_indel_model // string
    sentieon_dnascope_model           // channel

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
                intervals_name:intervals_.simpleName,
                variantcaller:'sentieon_dnascope'],
            cram_,
            crai,
            intervals_
            ]
        }

    emit_mode_items = sentieon_dnascope_emit_mode.split(',').each{ it -> it.toLowerCase().trim() }
    lst = emit_mode_items - 'gvcf'
    emit_vcf = lst.size() > 0 ? lst[0] : ''

    SENTIEON_DNASCOPE(
        cram_intervals_for_sentieon,
        fasta,
        fasta_fai,
        dbsnp.map{it -> [[:], it]},
        dbsnp_tbi.map{it -> [[:], it]},
        sentieon_dnascope_model.map{it -> [[:], it]},
        sentieon_dnascope_pcr_indel_model,
        emit_vcf,
        emit_mode_items.any{ it.equals('gvcf') })

    if (joint_germline) {
        genotype_intervals = SENTIEON_DNASCOPE.out.gvcf
            .join(SENTIEON_DNASCOPE.out.gvcf_tbi, failOnMismatch: true)
            .join(cram_intervals_for_sentieon, failOnMismatch: true)
            .map{ meta, gvcf_, tbi, cram_, crai, intervals_ -> [ meta, gvcf_, tbi, intervals_ ] }
    }

    // Figure out if using intervals or no_intervals.
    // Strip `intervals_name` here so the branch output meta is the same one we'll
    // group on later — fixes the prior `interval_name` typo which was a no-op.
    dnascope_vcf_branch = SENTIEON_DNASCOPE.out.vcf.map{
            meta, vcf_ -> [ meta - meta.subMap('intervals_name'), vcf_]
        }
        .branch{
            intervals:    it[0].num_intervals > 1
            no_intervals: it[0].num_intervals <= 1
        }

    dnascope_vcf_tbi_branch = SENTIEON_DNASCOPE.out.vcf_tbi.map{
            meta, vcf_tbi -> [ meta - meta.subMap('intervals_name'), vcf_tbi]
        }
        .branch{
            intervals:    it[0].num_intervals > 1
            no_intervals: it[0].num_intervals <= 1
        }

    haplotyper_gvcf_branch = SENTIEON_DNASCOPE.out.gvcf.map{
            meta, gvcf_ -> [ meta - meta.subMap('intervals_name'), gvcf_]
        }
        .branch{
            intervals:    it[0].num_intervals > 1
            no_intervals: it[0].num_intervals <= 1
        }

    haplotyper_gvcf_tbi_branch = SENTIEON_DNASCOPE.out.gvcf_tbi.map{
            meta, gvcf_tbi -> [ meta - meta.subMap('intervals_name'), gvcf_tbi]
        }
        .branch{
            intervals:    it[0].num_intervals > 1
            no_intervals: it[0].num_intervals <= 1
        }

    // Per-sample merge. Wrap the (already-`intervals_name`-stripped) meta in
    // groupKey *inside* the map so the GroupKey survives into groupTuple —
    // `groupKey` followed by a separate `meta - subMap(...)` would unwrap it
    // (Map.minus(Map) returns a plain LinkedHashMap), losing the size hint and
    // forcing groupTuple to wait until the upstream channel closes.
    vcfs_for_merging = dnascope_vcf_branch.intervals.map{
        meta, vcf_ -> [ groupKey(meta, meta.num_intervals), vcf_ ]
    }.groupTuple()

    // VCFs
    // Only when using intervals
    MERGE_SENTIEON_DNASCOPE_VCFS(vcfs_for_merging, dict)

    dnascope_vcf = Channel.empty().mix(
        MERGE_SENTIEON_DNASCOPE_VCFS.out.vcf,
        dnascope_vcf_branch.no_intervals)

    haplotyper_tbi = Channel.empty().mix(
        MERGE_SENTIEON_DNASCOPE_VCFS.out.tbi,
        dnascope_vcf_tbi_branch.no_intervals)

    // Remove no longer necessary field: num_intervals
    vcf = dnascope_vcf.map{ meta, vcf_ -> [ meta - meta.subMap('num_intervals'), vcf_ ] }
    vcf_tbi = haplotyper_tbi.map{ meta, tbi -> [ meta - meta.subMap('num_intervals'), tbi ] }

    // GVCFs — see vcfs_for_merging comment above for why the groupKey is built
    // inside the same map (preserves size hint into groupTuple).
    gvcfs_for_merging = haplotyper_gvcf_branch.intervals.map{
        meta, vcf_ -> [ groupKey(meta, meta.num_intervals), vcf_ ]
    }.groupTuple()

    MERGE_SENTIEON_DNASCOPE_GVCFS(gvcfs_for_merging, dict)

    gvcf = Channel.empty().mix(
        MERGE_SENTIEON_DNASCOPE_GVCFS.out.vcf,
        haplotyper_gvcf_branch.no_intervals)

    gvcf_tbi = Channel.empty().mix(
        MERGE_SENTIEON_DNASCOPE_GVCFS.out.tbi,
        haplotyper_gvcf_tbi_branch.no_intervals)

    versions = versions.mix(MERGE_SENTIEON_DNASCOPE_VCFS.out.versions)
    versions = versions.mix(MERGE_SENTIEON_DNASCOPE_GVCFS.out.versions)

    emit:
    versions
    vcf
    vcf_tbi
    gvcf
    gvcf_tbi
    genotype_intervals // For joint genotyping

}
