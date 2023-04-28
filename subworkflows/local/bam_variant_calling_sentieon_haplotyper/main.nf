include { BAM_MERGE_INDEX_SAMTOOLS                 } from '../bam_merge_index_samtools/main'
include { VCF_VARIANT_FILTERING_GATK               } from '../vcf_variant_filtering_gatk/main'
include { GATK4_HAPLOTYPECALLER                    } from '../../../modules/nf-core/gatk4/haplotypecaller/main'
include { GATK4_MERGEVCFS as MERGE_HAPLOTYPECALLER } from '../../../modules/nf-core/gatk4/mergevcfs/main'
include { SENTIEON_HAPLOTYPER                      } from '../../../modules/nf-core/sentieon/haplotyper/main'

include { BAM_JOINT_CALLING_GERMLINE_GATK                } from '../bam_joint_calling_germline_gatk/main'  // TO-DO: Figure out whether this is really needed!!!
include { GATK4_MERGEVCFS as MERGE_SENTIEON_HAPLOTYPER_VCFS  } from '../../../modules/nf-core/gatk4/mergevcfs/main'  // TO-DO: Figure out whether this is really needed!!!
include { GATK4_MERGEVCFS as MERGE_SENTIEON_HAPLOTYPER_GVCFS } from '../../../modules/nf-core/gatk4/mergevcfs/main'  // TO-DO: Figure out whether this is really needed!!!

workflow BAM_VARIANT_CALLING_SENTIEON_HAPLOTYPER {
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
    skip_haplotyper_filter       // boolean: [mandatory] [default: false] skip haplotyper filter

    main:
    versions = Channel.empty()

    gvcf = Channel.empty()  // SENTIEON

    vcf = Channel.empty()
    realigned_bam = Channel.empty()

    // Combine cram and intervals for spread and gather strategy
    cram_intervals = cram.combine(intervals)
        // Move num_intervals to meta map
        .map{ meta, cram, crai, intervals, num_intervals -> [ meta + [ num_intervals:num_intervals ], cram, crai, intervals, [] ] }

    // SENTIEON
    cram_intervals_for_sentieon = cram.combine(intervals)
        // Move num_intervals to meta map
        .map{ meta, cram, crai, intervals, num_intervals -> [ meta + [ num_intervals:num_intervals ], cram, crai, intervals ] }

    SENTIEON_HAPLOTYPER(
        'vcf',
        cram_intervals_for_sentieon,
        fasta,
        fasta_fai,
        dbsnp,
        dbsnp_tbi)

    versions = versions.mix(SENTIEON_HAPLOTYPER.out.versions)

/*
    // For joint genotyping // Not working for Sention!!! Err msg : "Join mismatch for the following entries.". Anyway, "genotype_intervals" may not be needed for Sentieon.
    genotype_intervals_for_sentieon = SENTIEON_HAPLOTYPER.out.gvcf
        .join(SENTIEON_HAPLOTYPER.out.gvcf_tbi, failOnMismatch: true)
        .join(cram_intervals_for_sentieon, failOnMismatch: true)
        .map{ meta, gvcf, tbi, cram, crai, intervals -> [ meta, gvcf, tbi, intervals ] }

    cram_intervals_for_sentieon.view{" cram_intervals_for_sentieon : $it"}

    SENTIEON_HAPLOTYPER.out.gvcf.view{ "cram_intervals_for_sentieon : $it"}

    genotype_intervals_for_sentieon = SENTIEON_HAPLOTYPER.out.gvcf
        .join(SENTIEON_HAPLOTYPER.out.gvcf_tbi, failOnMismatch: true)

    genotype_intervals_for_sentieon.view{" genotype_intervals_for_sentieon : $it"}
*/

    // Figure out if using intervals or no_intervals
    SENTIEON_HAPLOTYPER.out.vcf.branch{
            intervals:    it[0].num_intervals > 1
            no_intervals: it[0].num_intervals <= 1
        }.set{haplotypecaller_vcf_branch}

    SENTIEON_HAPLOTYPER.out.vcf_tbi.branch{
            intervals:    it[0].num_intervals > 1
            no_intervals: it[0].num_intervals <= 1
        }.set{haplotypecaller_vcf_tbi_branch}

    SENTIEON_HAPLOTYPER.out.gvcf.branch{
            intervals:    it[0].num_intervals > 1
            no_intervals: it[0].num_intervals <= 1
        }.set{haplotypecaller_gvcf_branch}

    SENTIEON_HAPLOTYPER.out.gvcf_tbi.branch{
            intervals:    it[0].num_intervals > 1
            no_intervals: it[0].num_intervals <= 1
        }.set{haplotypecaller_gvcf_tbi_branch}

    if (params.joint_germline) {
        // TO-DO:  Test and fix the joint_germline-workflow
        // merge vcf and tbis
        genotype_gvcf_to_call = Channel.empty().mix(SENTIEON_HAPLOTYPER.out.gvcf
                                                    .join(SENTIEON_HAPLOTYPER.out.gvcf_tbi)
                                                    .join(cram).map{ meta, gvcf, gvcf_tbi, cram, crai, intervals, dragstr_model ->
                                                            [ meta, gvcf, gvcf_tbi, intervals ]
                                                    })
        // make channels from labels
        dbsnp_vqsr        = params.dbsnp_vqsr        ? Channel.value(params.dbsnp_vqsr)        : Channel.empty()
        known_indels_vqsr = params.known_indels_vqsr ? Channel.value(params.known_indels_vqsr) : Channel.empty()
        known_snps_vqsr   = params.known_snps_vqsr   ? Channel.value(params.known_snps_vqsr)   : Channel.empty()


        BAM_JOINT_CALLING_GERMLINE_GATK(
            genotype_gvcf_to_call,
            fasta,
            fasta_fai,
            dict,
            dbsnp,
            dbsnp_tbi,
            dbsnp_vqsr,
            known_sites_indels,
            known_sites_indels_tbi,
            known_indels_vqsr,
            known_sites_snps,
            known_sites_snps_tbi,
            known_snps_vqsr)

        versions = versions.mix(BAM_JOINT_CALLING_GERMLINE_GATK.out.versions)

        gvcf = BAM_JOINT_CALLING_GERMLINE_GATK.out.genotype_vcf
    } else {

        // VCFs
        // Only when using intervals
        MERGE_SENTIEON_HAPLOTYPER_VCFS(
            haplotypecaller_vcf_branch.intervals
            .map{ meta, vcf ->

                new_meta = [
                                id:             meta.sample,
                                num_intervals:  meta.num_intervals,
                                patient:        meta.patient,
                                sample:         meta.sample,
                                sex:            meta.sex,
                                status:         meta.status
                            ]

                    [groupKey(new_meta, new_meta.num_intervals), vcf]
                }.groupTuple(),
            dict.map{it -> [ [ id:it.baseName ], it ]})


        versions = versions.mix(MERGE_SENTIEON_HAPLOTYPER_VCFS.out.versions)

        haplotypecaller_vcf = Channel.empty().mix(
            MERGE_SENTIEON_HAPLOTYPER_VCFS.out.vcf,
            haplotypecaller_vcf_branch.no_intervals)

        haplotypecaller_tbi = Channel.empty().mix(
            MERGE_SENTIEON_HAPLOTYPER_VCFS.out.tbi,
            haplotypecaller_vcf_tbi_branch.no_intervals)

        if (!skip_haplotyper_filter) {
            VCF_VARIANT_FILTERING_GATK(haplotypecaller_vcf.join(haplotypecaller_tbi),
                        fasta,
                        fasta_fai,
                        dict,
                        intervals_bed_combined,
                        known_sites_indels.concat(known_sites_snps).flatten().unique().collect(),
                        known_sites_indels_tbi.concat(known_sites_snps_tbi).flatten().unique().collect())

            vcf = VCF_VARIANT_FILTERING_GATK.out.filtered_vcf.map{ meta, vcf-> [[patient:meta.patient, sample:meta.sample, status:meta.status, sex:meta.sex, id:meta.sample, num_intervals:meta.num_intervals, variantcaller:"sentieon_haplotyper"], vcf]}
            versions = versions.mix(VCF_VARIANT_FILTERING_GATK.out.versions)

        } else vcf = haplotypecaller_vcf


        // add variantcaller to meta map and remove no longer necessary field: num_intervals
        vcf = vcf.map{ meta, vcf -> [ meta - meta.subMap('num_intervals') + [ variantcaller:'sentieon_haplotyper' ], vcf ] }


        // GVFs
        // Only when using intervals
        MERGE_SENTIEON_HAPLOTYPER_GVCFS(
            haplotypecaller_gvcf_branch.intervals
            .map{ meta, gvcf ->

                new_meta = [
                                id:             meta.sample,
                                num_intervals:  meta.num_intervals,
                                patient:        meta.patient,
                                sample:         meta.sample,
                                sex:            meta.sex,
                                status:         meta.status,
                                variantcaller:  "haplotyper"
                            ]

                    [groupKey(new_meta, new_meta.num_intervals), gvcf]
                }.groupTuple(),
            dict)

        versions = versions.mix(MERGE_SENTIEON_HAPLOTYPER_GVCFS.out.versions)

        gvcf = Channel.empty().mix(
            MERGE_SENTIEON_HAPLOTYPER_GVCFS.out.vcf,
            haplotypecaller_gvcf_branch.no_intervals)


    }



/*
    // The original GATK-haplotypecaller code. Just here for comparison during development.

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

    /* if (!skip_haplotypecaller_filter) {     // Temporary disabled */
/*
    if (false) {

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
*/
    emit:
    versions
    vcf
    gvcf

}
