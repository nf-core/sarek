//
// GATK4 HAPLOTYPACALLER germline variant calling:
//
// For all modules here:
// A when clause condition is defined in the conf/modules.config to determine if the module should be run

include { BAM_MERGE_INDEX_SAMTOOLS                            } from '../bam_merge_index_samtools/main'
include { GATK4_HAPLOTYPECALLER                               } from '../../../modules/nf-core/gatk4/haplotypecaller/main'
include { GATK4_MERGEVCFS            as MERGE_HAPLOTYPECALLER } from '../../../modules/nf-core/gatk4/mergevcfs/main'

workflow BAM_VARIANT_CALLING_HAPLOTYPECALLER {
    take:
    cram                         // channel: [mandatory] [ meta, cram, crai, interval.bed ]
    fasta                        // channel: [mandatory]
    fasta_fai                    // channel: [mandatory]
    dict                         // channel: [mandatory]
    dbsnp                        // channel: [optional]
    dbsnp_tbi                    // channel: [optional]
    intervals                    // channel: [mandatory] [ intervals, num_intervals ] or [ [], 0 ] if no intervals

    main:
    versions = Channel.empty()

    vcf           = Channel.empty()
    realigned_bam = Channel.empty()

    // Combine cram and intervals for spread and gather strategy
    cram_intervals = cram.combine(intervals)
        // Move num_intervals to meta map
        // Add interval_name to allow correct merging with interval files
        .map{ meta, cram, crai, intervals, num_intervals -> [ meta + [ interval_name:intervals.baseName, num_intervals:num_intervals, variantcaller:'haplotypecaller' ], cram, crai, intervals, [] ] }

    GATK4_HAPLOTYPECALLER(
        cram_intervals,
        fasta.map{ it -> [[:], it] },
        fasta_fai.map{ it -> [[:], it] },
        dict,
        dbsnp.map{ it -> [[:], it] },
        dbsnp_tbi.map{ it -> [[:], it] })

    // For joint genotyping
    gvcf_tbi_intervals = GATK4_HAPLOTYPECALLER.out.vcf
        .join(GATK4_HAPLOTYPECALLER.out.tbi, failOnMismatch: true)
        .join(cram_intervals, failOnMismatch: true)
        .map{ meta, gvcf, tbi, cram, crai, intervals, dragstr_model -> [ meta, gvcf, tbi, intervals ] }

    // Figuring out if there is one or more vcf(s) from the same sample
    haplotypecaller_vcf = GATK4_HAPLOTYPECALLER.out.vcf.map{
            meta, vcf -> [ meta - meta.subMap('interval_name'), vcf]
        }
        .branch{
        // Use meta.num_intervals to asses number of intervals
            intervals:    it[0].num_intervals > 1
            no_intervals: it[0].num_intervals <= 1
        }

    // Figuring out if there is one or more tbi(s) from the same sample
    haplotypecaller_tbi = GATK4_HAPLOTYPECALLER.out.tbi.map{
            meta, tbi -> [ meta - meta.subMap('interval_name'), tbi]
        }.branch{
        // Use meta.num_intervals to asses number of intervals
            intervals:    it[0].num_intervals > 1
            no_intervals: it[0].num_intervals <= 1
        }

    // Figuring out if there is one or more bam(s) from the same sample
    haplotypecaller_bam = GATK4_HAPLOTYPECALLER.out.bam.map{
            meta, bam -> [ meta - meta.subMap('interval_name'), bam]
        }.branch{
        // Use meta.num_intervals to asses number of intervals
            intervals:    it[0].num_intervals > 1
            no_intervals: it[0].num_intervals <= 1
        }

    // Only when using intervals
    MERGE_HAPLOTYPECALLER(haplotypecaller_vcf.intervals.map{ meta, vcf -> [ groupKey(meta, meta.num_intervals), vcf ] }.groupTuple(), dict)

    haplotypecaller_vcf = Channel.empty().mix(
            MERGE_HAPLOTYPECALLER.out.vcf,
            haplotypecaller_vcf.no_intervals)

    haplotypecaller_tbi = Channel.empty().mix(
            MERGE_HAPLOTYPECALLER.out.tbi,
            haplotypecaller_tbi.no_intervals)

    // BAM output
    BAM_MERGE_INDEX_SAMTOOLS(haplotypecaller_bam.intervals
        .map{ meta, bam -> [ groupKey(meta, meta.num_intervals), bam ] }
        .groupTuple()
        .mix(haplotypecaller_bam.no_intervals))

    realigned_bam = BAM_MERGE_INDEX_SAMTOOLS.out.bam_bai

    versions = versions.mix(GATK4_HAPLOTYPECALLER.out.versions)
    versions = versions.mix(MERGE_HAPLOTYPECALLER.out.versions)

    // Remove no longer necessary field: num_intervals
    vcf = haplotypecaller_vcf.map{ meta, vcf -> [ meta - meta.subMap('num_intervals'), vcf ] }
    tbi = haplotypecaller_tbi.map{ meta, tbi -> [ meta - meta.subMap('num_intervals'), tbi ] }

    emit:
    gvcf_tbi_intervals // For joint genotyping
    realigned_bam      // Optional
    vcf                // vcf
    tbi                // tbi

    versions
}
