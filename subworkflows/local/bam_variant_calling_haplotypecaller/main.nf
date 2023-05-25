include { BAM_MERGE_INDEX_SAMTOOLS                 } from '../bam_merge_index_samtools/main'
include { VCF_VARIANT_FILTERING_GATK               } from '../vcf_variant_filtering_gatk/main'
include { GATK4_HAPLOTYPECALLER                    } from '../../../modules/nf-core/gatk4/haplotypecaller/main'
include { GATK4_MERGEVCFS as MERGE_HAPLOTYPECALLER } from '../../../modules/nf-core/gatk4/mergevcfs/main'

workflow BAM_VARIANT_CALLING_HAPLOTYPECALLER {
    take:
    cram                         // channel: [mandatory] [ meta, cram, crai, interval.bed ]
    fasta                        // channel: [mandatory]
    fasta_fai                    // channel: [mandatory]
    dict                         // channel: [mandatory]
    dbsnp                        // channel: [optional]
    dbsnp_tbi                    // channel: [optional]
    dbsnp_vqsr                   // channel: [optional]
    known_sites_indels           // channel: [optional]
    known_sites_indels_tbi       // channel: [optional]
    known_indels_vqsr            // channel: [optional]
    known_sites_snps             // channel: [optional]
    known_sites_snps_tbi         // channel: [optional]
    known_snps_vqsr              // channel: [optional]
    intervals                    // channel: [mandatory] [ intervals, num_intervals ] or [ [], 0 ] if no intervals
    intervals_bed_combined       // channel: [mandatory] intervals/target regions in one file unzipped, no_intervals.bed if no_intervals
    skip_haplotypecaller_filter  // boolean: [mandatory] [default: false] skip haplotypecaller filter

    main:
    versions = Channel.empty()

    vcf = Channel.empty()
    realigned_bam = Channel.empty()

    // Combine cram and intervals for spread and gather strategy
    cram_intervals = cram.combine(intervals)
        // Move num_intervals to meta map
        .map{ meta, cram, crai, intervals, num_intervals -> [ meta + [ num_intervals:num_intervals ], cram, crai, intervals, [] ] }


    GATK4_HAPLOTYPECALLER(cram_intervals, fasta, fasta_fai, dict, dbsnp, dbsnp_tbi)

    // For joint genotyping
    genotype_intervals = GATK4_HAPLOTYPECALLER.out.vcf
        .join(GATK4_HAPLOTYPECALLER.out.tbi, failOnMismatch: true)
        .join(cram_intervals, failOnMismatch: true)
        .map{ meta, gvcf, tbi, cram, crai, intervals, dragstr_model -> [ meta, gvcf, tbi, intervals ] }

    // Figuring out if there is one or more vcf(s) from the same sample
    haplotypecaller_vcf = GATK4_HAPLOTYPECALLER.out.vcf.branch{
        // Use meta.num_intervals to asses number of intervals
            intervals:    it[0].num_intervals > 1
            no_intervals: it[0].num_intervals <= 1
        }

    // Figuring out if there is one or more tbi(s) from the same sample
    haplotypecaller_tbi = GATK4_HAPLOTYPECALLER.out.tbi.branch{
        // Use meta.num_intervals to asses number of intervals
            intervals:    it[0].num_intervals > 1
            no_intervals: it[0].num_intervals <= 1
        }

    // Figuring out if there is one or more bam(s) from the same sample
    haplotypecaller_bam = GATK4_HAPLOTYPECALLER.out.bam.branch{
        // Use meta.num_intervals to asses number of intervals
            intervals:    it[0].num_intervals > 1
            no_intervals: it[0].num_intervals <= 1
        }

    // Only when using intervals
    MERGE_HAPLOTYPECALLER(haplotypecaller_vcf.intervals
        .map{ meta, vcf -> [ groupKey(meta, meta.num_intervals), vcf ] }
        .groupTuple(),
        dict.map{ it -> [ [ id:'dict' ], it ] })

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

    if (!skip_haplotypecaller_filter) {

        VCF_VARIANT_FILTERING_GATK(
            haplotypecaller_vcf.join(haplotypecaller_tbi, failOnDuplicate: true, failOnMismatch: true),
            fasta,
            fasta_fai,
            dict,
            intervals_bed_combined,
            known_sites_indels.concat(known_sites_snps).flatten().unique().collect(),
            known_sites_indels_tbi.concat(known_sites_snps_tbi).flatten().unique().collect())

        vcf = VCF_VARIANT_FILTERING_GATK.out.filtered_vcf

        versions = versions.mix(VCF_VARIANT_FILTERING_GATK.out.versions)

    } else vcf = haplotypecaller_vcf

    versions = versions.mix(GATK4_HAPLOTYPECALLER.out.versions)
    versions = versions.mix(MERGE_HAPLOTYPECALLER.out.versions)

    // add variantcaller to meta map and remove no longer necessary field: num_intervals
    vcf = vcf.map{ meta, vcf -> [ meta - meta.subMap('num_intervals') + [ variantcaller:'haplotypecaller' ], vcf ] }

    emit:
    genotype_intervals // For joint genotyping
    realigned_bam      // Optionnal
    vcf                // vcf filtered or not

    versions
}
