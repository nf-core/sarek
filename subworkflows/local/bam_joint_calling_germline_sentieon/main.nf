//
// JOINT GERMLINE CALLING
//
// Merge samples perform joint genotyping with SENTIEON_GVCFTYPER
//

include { BCFTOOLS_SORT                                      } from '../../../modules/nf-core/bcftools/sort/main'
include { GATK4_MERGEVCFS      as MERGE_GENOTYPEGVCFS        } from '../../../modules/nf-core/gatk4/mergevcfs/main'
include { GATK4_MERGEVCFS      as MERGE_VQSR                 } from '../../../modules/nf-core/gatk4/mergevcfs/main'
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

    main:
    versions = Channel.empty()

    sentieon_input = input
        .map{ meta, gvcf, tbi, intervals -> [ [ id:'joint_variant_calling', intervals_name:intervals.simpleName, num_intervals:meta.num_intervals ], gvcf, tbi, intervals ] }
        .groupTuple(by:[0, 3])

    SENTIEON_GVCFTYPER(sentieon_input, fasta, fai, dbsnp, dbsnp_tbi)

    BCFTOOLS_SORT(SENTIEON_GVCFTYPER.out.vcf_gz)

    gvcf_to_merge = BCFTOOLS_SORT.out.vcf.map{ meta, vcf -> [ meta.subMap('num_intervals') + [ id:'joint_variant_calling', patient:'all_samples', variantcaller:'sentieon_haplotyper' ], vcf ]}.groupTuple()

    // Merge scatter/gather vcfs & index
    // Rework meta for variantscalled.csv and annotation tools
    MERGE_GENOTYPEGVCFS(gvcf_to_merge, dict)

    vqsr_input = MERGE_GENOTYPEGVCFS.out.vcf.join(MERGE_GENOTYPEGVCFS.out.tbi, failOnDuplicate: true)
    indels_resource_label = known_indels_vqsr.mix(dbsnp_vqsr).collect()
    snps_resource_label = known_snps_vqsr.mix(dbsnp_vqsr).collect()

    // Recalibrate INDELs and SNPs separately
    SENTIEON_VARCAL_INDEL(
        vqsr_input,
        resource_indels_vcf,
        resource_indels_tbi,
        indels_resource_label,
        fasta,
        fai)

    SENTIEON_VARCAL_SNP(
        vqsr_input,
        resource_snps_vcf,
        resource_snps_tbi,
        snps_resource_label,
        fasta,
        fai)

    //Prepare INDELs and SNPs separately for Sentieons applyvarcal

    // Join results of variant recalibration into a single channel tuple
    // Rework meta for variantscalled.csv and annotation tools
    vqsr_input_snp = vqsr_input.join(SENTIEON_VARCAL_SNP.out.recal, failOnDuplicate: true)
        .join(SENTIEON_VARCAL_SNP.out.idx, failOnDuplicate: true)
        .join(SENTIEON_VARCAL_SNP.out.tranches, failOnDuplicate: true)
        .map{ meta, vcf, tbi, recal, index, tranche -> [ meta + [ id:'recalibrated_joint_variant_calling' ], vcf, tbi, recal, index, tranche ] }

    // Join results of variant recalibration into a single channel tuple
    // Rework meta for variantscalled.csv and annotation tools
    vqsr_input_indel = vqsr_input.join(SENTIEON_VARCAL_INDEL.out.recal, failOnDuplicate: true)
        .join(SENTIEON_VARCAL_INDEL.out.idx, failOnDuplicate: true)
        .join(SENTIEON_VARCAL_INDEL.out.tranches, failOnDuplicate: true)
        .map{ meta, vcf, tbi, recal, index, tranche -> [ meta + [ id:'recalibrated_joint_variant_calling' ], vcf, tbi, recal, index, tranche ] }

    SENTIEON_APPLYVARCAL_SNP(
        vqsr_input_snp,
        fasta,
        fai)

    SENTIEON_APPLYVARCAL_INDEL(
        vqsr_input_indel,
        fasta,
        fai)

    vqsr_snp_vcf = SENTIEON_APPLYVARCAL_SNP.out.vcf
    vqsr_indel_vcf = SENTIEON_APPLYVARCAL_INDEL.out.vcf

    //Merge VQSR outputs into final VCF
    MERGE_VQSR(
        vqsr_snp_vcf.mix(vqsr_indel_vcf).groupTuple(),
        dict
    )

    genotype_vcf   = MERGE_GENOTYPEGVCFS.out.vcf.mix(MERGE_VQSR.out.vcf)
    genotype_index = MERGE_GENOTYPEGVCFS.out.tbi.mix(MERGE_VQSR.out.tbi)

    versions = versions.mix(SENTIEON_GVCFTYPER.out.versions)
    versions = versions.mix(SENTIEON_VARCAL_SNP.out.versions)
    versions = versions.mix(SENTIEON_VARCAL_INDEL.out.versions)
    versions = versions.mix(SENTIEON_APPLYVARCAL_INDEL.out.versions)

    emit:
    genotype_index  // channel: [ val(meta), [ tbi ] ]
    genotype_vcf    // channel: [ val(meta), [ vcf ] ]

    versions        // channel: [ versions.yml ]
}
