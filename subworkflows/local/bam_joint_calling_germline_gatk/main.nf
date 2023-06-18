//
// merge samples with genomicsdbimport, perform joint genotyping with genotypeGVCFS
include { BCFTOOLS_SORT                                          } from '../../../modules/nf-core/bcftools/sort/main'
include { GATK4_APPLYVQSR as GATK4_APPLYVQSR_INDEL               } from '../../../modules/nf-core/gatk4/applyvqsr/main'
include { GATK4_APPLYVQSR as GATK4_APPLYVQSR_SNP                 } from '../../../modules/nf-core/gatk4/applyvqsr/main'
include { GATK4_GENOMICSDBIMPORT                                 } from '../../../modules/nf-core/gatk4/genomicsdbimport/main'
include { GATK4_GENOTYPEGVCFS                                    } from '../../../modules/nf-core/gatk4/genotypegvcfs/main'
include { GATK4_MERGEVCFS as MERGE_GENOTYPEGVCFS                 } from '../../../modules/nf-core/gatk4/mergevcfs/main'
include { GATK4_MERGEVCFS as MERGE_VQSR                          } from '../../../modules/nf-core/gatk4/mergevcfs/main'
include { GATK4_VARIANTRECALIBRATOR as VARIANTRECALIBRATOR_INDEL } from '../../../modules/nf-core/gatk4/variantrecalibrator/main'
include { GATK4_VARIANTRECALIBRATOR as VARIANTRECALIBRATOR_SNP   } from '../../../modules/nf-core/gatk4/variantrecalibrator/main'

workflow BAM_JOINT_CALLING_GERMLINE_GATK {
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
    genomicsdb_workspace

    main:
    versions = Channel.empty()

    // Map input for GenomicsDBImport
    // Rename based on num_intervals, group all samples by their interval_name/interval_file and restructure for channel
    // Group by [0, 3] to avoid a list of metas and make sure that any intervals

    gendb_input = input
        .map{ meta, gvcf, tbi, intervals -> [ [ id:'joint_variant_calling', intervals_name:intervals.simpleName, num_intervals:meta.num_intervals ], gvcf, tbi, intervals ] }
        .groupTuple(by:[0, 3])
        .map{ meta, gvcf, tbi, intervals -> [ meta, gvcf, tbi, intervals, [] ] }
        .combine(genomicsdb_workspace)

    // Convert all sample vcfs into a genomicsdb workspace using genomicsdbimport
    // If genomicsdb_workspace, then the boolean will evaluate to true, hacky work-around for the module
    GATK4_GENOMICSDBIMPORT(gendb_input, false, genomicsdb_workspace, false)

    genotype_input = GATK4_GENOMICSDBIMPORT.out.genomicsdb.map{ meta, genomicsdb -> [ meta, genomicsdb, [], [], [] ] }

    // Joint genotyping performed using GenotypeGVCFs
    // Sort vcfs called by interval within each VCF

    GATK4_GENOTYPEGVCFS(genotype_input, fasta, fai, dict.map{ meta, dict -> [ dict ] }, dbsnp, dbsnp_tbi)

    BCFTOOLS_SORT(GATK4_GENOTYPEGVCFS.out.vcf)
    gvcf_to_merge = BCFTOOLS_SORT.out.vcf.map{ meta, vcf -> [ meta.subMap('num_intervals') + [ id:'joint_variant_calling', patient:'all_samples', variantcaller:'haplotypecaller' ], vcf ]}.groupTuple()

    // Merge scatter/gather vcfs & index
    // Rework meta for variantscalled.csv and annotation tools
    MERGE_GENOTYPEGVCFS(gvcf_to_merge, dict)

    vqsr_input = MERGE_GENOTYPEGVCFS.out.vcf.join(MERGE_GENOTYPEGVCFS.out.tbi, failOnDuplicate: true)
    indels_resource_label = known_indels_vqsr.mix(dbsnp_vqsr).collect()
    snps_resource_label = known_snps_vqsr.mix(dbsnp_vqsr).collect()

    // Recalibrate INDELs and SNPs separately
    VARIANTRECALIBRATOR_INDEL(
        vqsr_input,
        resource_indels_vcf,
        resource_indels_tbi,
        indels_resource_label,
        fasta,
        fai,
        dict.map{ meta, dict -> [ dict ] })

    VARIANTRECALIBRATOR_SNP(
        vqsr_input,
        resource_snps_vcf,
        resource_snps_tbi,
        snps_resource_label,
        fasta,
        fai,
        dict.map{ meta, dict -> [ dict ] })

    //Prepare INDELs and SNPs separately for ApplyVQSR

    // Join results of variant recalibration into a single channel tuple
    // Rework meta for variantscalled.csv and annotation tools
    vqsr_input_snp = vqsr_input.join(VARIANTRECALIBRATOR_SNP.out.recal, failOnDuplicate: true)
        .join(VARIANTRECALIBRATOR_SNP.out.idx, failOnDuplicate: true)
        .join(VARIANTRECALIBRATOR_SNP.out.tranches, failOnDuplicate: true)
        .map{ meta, vcf, tbi, recal, index, tranche -> [ meta - meta.subMap('id') + [ id:'recalibrated_joint_variant_calling' ], vcf, tbi, recal, index, tranche ] }

    // Join results of variant recalibration into a single channel tuple
    // Rework meta for variantscalled.csv and annotation tools
    vqsr_input_indel = vqsr_input.join(VARIANTRECALIBRATOR_INDEL.out.recal, failOnDuplicate: true)
        .join(VARIANTRECALIBRATOR_INDEL.out.idx, failOnDuplicate: true)
        .join(VARIANTRECALIBRATOR_INDEL.out.tranches, failOnDuplicate: true)
        .map{ meta, vcf, tbi, recal, index, tranche -> [ meta - meta.subMap('id') + [ id:'recalibrated_joint_variant_calling' ], vcf, tbi, recal, index, tranche ] }

    GATK4_APPLYVQSR_SNP(
        vqsr_input_snp,
        fasta,
        fai,
        dict.map{ meta, dict -> [ dict ] })

    GATK4_APPLYVQSR_INDEL(
        vqsr_input_indel,
        fasta,
        fai,
        dict.map{ meta, dict -> [ dict ] })

    vqsr_snp_vcf = GATK4_APPLYVQSR_SNP.out.vcf
    vqsr_indel_vcf = GATK4_APPLYVQSR_INDEL.out.vcf

    //Merge VQSR outputs into final VCF
    MERGE_VQSR(vqsr_snp_vcf.mix(vqsr_indel_vcf).groupTuple(), dict)

    genotype_vcf   = MERGE_GENOTYPEGVCFS.out.vcf.mix(MERGE_VQSR.out.vcf)
    genotype_index = MERGE_GENOTYPEGVCFS.out.tbi.mix(MERGE_VQSR.out.tbi)

    versions = versions.mix(GATK4_GENOMICSDBIMPORT.out.versions)
    versions = versions.mix(GATK4_GENOTYPEGVCFS.out.versions)
    versions = versions.mix(VARIANTRECALIBRATOR_SNP.out.versions)
    versions = versions.mix(GATK4_APPLYVQSR_SNP.out.versions)

    emit:
    genotype_index  // channel: [ val(meta), [ tbi ] ]
    genotype_vcf    // channel: [ val(meta), [ vcf ] ]

    versions        // channel: [ versions.yml ]
}
