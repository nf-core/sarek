//
// JOINT GERMLINE CALLING
//
// Merge samples perform joint genotyping with SENTIEON_GVCFTYPER
//

include { BCFTOOLS_SORT                                      } from '../../../modules/nf-core/bcftools/sort/main'
include { GATK4_MERGEVCFS      as MERGE_GENOTYPEGVCFS        } from '../../../modules/nf-core/gatk4/mergevcfs/main'
include { SENTIEON_APPLYVARCAL as SENTIEON_APPLYVARCAL_INDEL } from '../../../modules/nf-core/sentieon/applyvarcal/main'
include { SENTIEON_APPLYVARCAL as SENTIEON_APPLYVARCAL_SNP   } from '../../../modules/nf-core/sentieon/applyvarcal/main'
include { SENTIEON_GVCFTYPER                                 } from '../../../modules/nf-core/sentieon/gvcftyper/main'
include { SENTIEON_VARCAL      as SENTIEON_VARCAL_INDEL      } from '../../../modules/nf-core/sentieon/varcal/main'
include { SENTIEON_VARCAL      as SENTIEON_VARCAL_SNP        } from '../../../modules/nf-core/sentieon/varcal/main'

workflow BAM_JOINT_CALLING_GERMLINE_SENTIEON {
    take:
    input                // channel: [ meta, [ input ], [ input_index ], intervals ]
    fasta                // channel: [ fasta ]
    fai                  // channel: [ fasta_fai ]
    dict                 // channel: [ dict ]
    dbsnp
    dbsnp_tbi
    dbsnp_vqsr
    resource_indels_vcf
    resource_indels_tbi
    known_indels_vqsr
    resource_snps_vcf
    resource_snps_tbi
    known_snps_vqsr
    variant_caller

    main:
    versions = Channel.empty()

    sentieon_input = input
        .map{ meta, gvcf, tbi, intervals -> [ [ id:'joint_variant_calling', intervals_name:intervals.baseName, num_intervals:meta.num_intervals ], gvcf, tbi, intervals ] }
        .groupTuple(by:[0, 3])

    SENTIEON_GVCFTYPER(sentieon_input, fasta, fai, dbsnp, dbsnp_tbi)

    BCFTOOLS_SORT(SENTIEON_GVCFTYPER.out.vcf_gz)

    gvcf_to_merge = BCFTOOLS_SORT.out.vcf.map{ meta, vcf -> [ meta.subMap('num_intervals') + [ id:'joint_variant_calling', patient:'all_samples', variantcaller:variant_caller ], vcf ]}.groupTuple()

    // Merge scatter/gather vcfs & index
    // Rework meta for variantscalled.csv and annotation tools
    MERGE_GENOTYPEGVCFS(gvcf_to_merge, dict)

    merged_vcf = MERGE_GENOTYPEGVCFS.out.vcf.map{meta, vcf -> [[id: 'joint_variant_calling'] , vcf]}
    merged_tbi = MERGE_GENOTYPEGVCFS.out.tbi.map{meta, tbi -> [[id: 'joint_variant_calling'] , tbi]}

    if (variant_caller == 'sentieon_dnascope') {
        // As advised by Don Freed (Sentieon), VQSR is skipped for DnaScope
        genotype_vcf   = merged_vcf.map{
            meta, vcf -> [ meta + [ patient:"all_samples", variantcaller:'sentieon_dnascope'], vcf ]
        }
        genotype_index = merged_tbi.map{
            meta, tbi -> [ meta + [ patient:"all_samples", variantcaller:'sentieon_dnascope'], tbi ]
        }
    } else {
        vqsr_input = MERGE_GENOTYPEGVCFS.out.vcf.join(MERGE_GENOTYPEGVCFS.out.tbi, failOnDuplicate: true)
        indels_resource_label = known_indels_vqsr.mix(dbsnp_vqsr).collect()
        snps_resource_label = known_snps_vqsr.mix(dbsnp_vqsr).collect()

        // Recalibrate INDELs and SNPs separately
        SENTIEON_VARCAL_INDEL(
            vqsr_input,
            resource_indels_vcf,
            resource_indels_tbi,
            indels_resource_label,
            fasta.map{meta, it -> [ it ]},
            fai.map{meta, it -> [ it ]})

        SENTIEON_VARCAL_SNP(
            vqsr_input,
            resource_snps_vcf,
            resource_snps_tbi,
            snps_resource_label,
            fasta.map{meta, it -> [ it ]},
            fai.map{meta, it -> [ it ]})

        //Prepare SNPs and INDELs for Sentieon's applyvarcal
        // Step 1. : applyvarcal to SNPs
        // Step 2. : Use SENTIEON_APPLYVARCAL_SNP output and run ApplyVQSR_INDEL. This avoids duplicate entries in the vcf as described here: https://hpc.nih.gov/training/gatk_tutorial/vqsr.html

        // Join results of variant recalibration into a single channel tuple
        // Rework meta for variantscalled.csv and annotation tools
        vqsr_input_snp = vqsr_input.join(SENTIEON_VARCAL_SNP.out.recal, failOnDuplicate: true)
            .join(SENTIEON_VARCAL_SNP.out.idx, failOnDuplicate: true)
            .join(SENTIEON_VARCAL_SNP.out.tranches, failOnDuplicate: true)
            .map{ meta, vcf, tbi, recal, index, tranche -> [ meta + [ id:'recalibrated_joint_variant_calling' ], vcf, tbi, recal, index, tranche ] }

        SENTIEON_APPLYVARCAL_SNP(vqsr_input_snp, fasta, fai)

        // Join results of SENTIEON_APPLYVARCAL_SNP and use as input for SENTIEON_APPLYVARCAL_INDEL to avoid duplicate entries in the result
        // Rework meta for variantscalled.csv and annotation tools
        vqsr_input_indel = SENTIEON_APPLYVARCAL_SNP.out.vcf.join(SENTIEON_APPLYVARCAL_SNP.out.tbi).map{ meta, vcf, tbi -> [ meta + [ id:'joint_variant_calling' ], vcf, tbi ]}
            .join(SENTIEON_VARCAL_INDEL.out.recal, failOnDuplicate: true)
            .join(SENTIEON_VARCAL_INDEL.out.idx, failOnDuplicate: true)
            .join(SENTIEON_VARCAL_INDEL.out.tranches, failOnDuplicate: true)
            .map{ meta, vcf, tbi, recal, index, tranche -> [ meta + [ id:'recalibrated_joint_variant_calling' ], vcf, tbi, recal, index, tranche ] }

        SENTIEON_APPLYVARCAL_INDEL(vqsr_input_indel, fasta, fai)

        // The following is an ugly monster to achieve the following:
        // When MERGE_GENOTYPEGVCFS and SENTIEON_APPLYVARCAL are run, then use output from SENTIEON_APPLYVARCAL
        // When MERGE_GENOTYPEGVCFS and NOT SENTIEON_APPLYVARCAL, then use the output from MERGE_GENOTYPEGVCFS

        // Remap for both to have the same key, if ApplyBQSR is not run, the channel is empty --> populate with empty elements
        vqsr_vcf_for_join = SENTIEON_APPLYVARCAL_INDEL.out.vcf.ifEmpty([[:], []]).map{meta, vcf -> [[id: 'joint_variant_calling'] , vcf]}
        vqsr_tbi_for_join = SENTIEON_APPLYVARCAL_INDEL.out.tbi.ifEmpty([[:], []]).map{meta, tbi -> [[id: 'joint_variant_calling'] , tbi]}

        // Join on metamap
        // If both --> meta, vcf_merged, vcf_bqsr
        // If not VQSR --> meta, vcf_merged, []
        // if the second is empty, use the first
        genotype_vcf = merged_vcf.join(vqsr_vcf_for_join, remainder: true).map{
            meta, joint_vcf, recal_vcf ->

            vcf_out = recal_vcf ?: joint_vcf

            [[id:"joint_variant_calling", patient:"all_samples", variantcaller:"sentieon_haplotyper"], vcf_out]
        }

        genotype_index = merged_tbi.join(vqsr_tbi_for_join, remainder: true).map{
            meta, joint_tbi, recal_tbi ->

            tbi_out = recal_tbi ?: joint_tbi
            [[id:"joint_variant_calling", patient:"all_samples", variantcaller:"sentieon_haplotyper"], tbi_out]
        }

        versions = versions.mix(SENTIEON_VARCAL_SNP.out.versions)
        versions = versions.mix(SENTIEON_VARCAL_INDEL.out.versions)
        versions = versions.mix(SENTIEON_APPLYVARCAL_INDEL.out.versions)
    }

    versions = versions.mix(SENTIEON_GVCFTYPER.out.versions)

    emit:
    genotype_index  // channel: [ val(meta), [ tbi ] ]
    genotype_vcf    // channel: [ val(meta), [ vcf ] ]

    versions        // channel: [ versions.yml ]
}
