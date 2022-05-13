//
// Recalibrate with variantrecalibrator & applyvqsr.
//

include { GATK4_VARIANTRECALIBRATOR as GATK4_VARIANTRECALIBRATOR_SNP } from '../../../modules/gatk4/variantrecalibrator/main'
include { GATK4_VARIANTRECALIBRATOR as GATK4_VARIANTRECALIBRATOR_INDEL } from '../../../modules/gatk4/variantrecalibrator/main'
include { GATK4_APPLYVQSR as GATK4_APPLYVQSR_SNP } from '../../../modules/gatk4/applyvqsr/main'
include { GATK4_APPLYVQSR as GATK4_APPLYVQSR_INDEL } from '../../../modules/gatk4/applyvqsr/main'
include { GATK4_SELECTVARIANTS as GATK4_SELECTVARIANTS_SNP } from '../../../modules/gatk4/selectvariants/main'
include { GATK4_SELECTVARIANTS as GATK4_SELECTVARIANTS_INDEL } from '../../../modules/gatk4/selectvariants/main'
include { GATK4_SELECTVARIANTS as GATK4_SELECTVARIANTS_NORECAL } from '../../../modules/gatk4/selectvariants/main'
include { GATK4_MERGEVCFS as GATK4_MERGEVCFS_RECALIBRATED } from '../../../modules/gatk4/mergevcfs/main'

workflow GATK_VQSR {
    take:
    input            // channel: [ val(meta), [ input ], [ input_index ], [] ]
    fasta            // channel: /path/to/reference/fasta
    fai              // channel: /path/to/reference/fasta/index
    dict             // channel: /path/to/reference/fasta/dictionary
    allelespecific   // channel: true/false run allelespecific mode of vqsr modules
    resources_SNP        // channel: [[resource, vcfs, forvariantrecal], [resource, tbis, forvariantrecal], [resource, labels, forvariantrecal]]
    resources_INDEL        // channel: [[resource, vcfs, forvariantrecal], [resource, tbis, forvariantrecal], [resource, labels, forvariantrecal]]
    annotation_SNP       // channel: [annotations, to, use, for, variantrecal, filtering]
    annotation_INDEL       // channel: [annotations, to, use, for, variantrecal, filtering]
    create_rscript   // channel: true/false whether to generate rscript plots in variantrecal
    truthsensitivity // channel: 0-100.0 truthsensitivity cutoff for applyvqsr

    main:
    ch_versions = Channel.empty()
    ch_select_variants_in = input

    GATK4_SELECTVARIANTS_SNP (ch_select_variants_in)
    ch_versions   = ch_versions.mix(GATK4_SELECTVARIANTS_SNP.out.versions)

    ch_vrecal_snp_in = GATK4_SELECTVARIANTS_SNP.out.vcf.combine(GATK4_SELECTVARIANTS_SNP.out.tbi, by: 0)
    GATK4_VARIANTRECALIBRATOR_SNP ( ch_vrecal_snp_in, fasta, fai, dict, allelespecific, resources_SNP, annotation_SNP, 'SNP', create_rscript )
    ch_versions   = ch_versions.mix(GATK4_VARIANTRECALIBRATOR_SNP.out.versions)

    ch_snp_recal      = GATK4_VARIANTRECALIBRATOR_SNP.out.recal
    ch_snp_idx        = GATK4_VARIANTRECALIBRATOR_SNP.out.idx
    ch_snp_tranches   = GATK4_VARIANTRECALIBRATOR_SNP.out.tranches
    ch_snp_vqsr_in    = ch_vrecal_snp_in.combine(ch_snp_recal, by: 0).combine(ch_snp_idx, by: 0).combine(ch_snp_tranches, by: 0)
    GATK4_APPLYVQSR_SNP ( ch_snp_vqsr_in, fasta, fai, dict, allelespecific, truthsensitivity, 'SNP' )
    ch_versions   = ch_versions.mix(GATK4_APPLYVQSR_SNP.out.versions)

    GATK4_SELECTVARIANTS_INDEL (ch_select_variants_in)

    ch_vrecal_indel_in = GATK4_SELECTVARIANTS_INDEL.out.vcf.combine(GATK4_SELECTVARIANTS_INDEL.out.tbi, by: 0)
    GATK4_VARIANTRECALIBRATOR_INDEL ( ch_vrecal_indel_in, fasta, fai, dict, allelespecific, resources_INDEL, annotation_INDEL, 'INDEL', create_rscript )

    ch_indel_recal      = GATK4_VARIANTRECALIBRATOR_INDEL.out.recal
    ch_indel_idx        = GATK4_VARIANTRECALIBRATOR_INDEL.out.idx
    ch_indel_tranches   = GATK4_VARIANTRECALIBRATOR_INDEL.out.tranches
    ch_indel_vqsr_in    = ch_vrecal_indel_in.combine(ch_indel_recal, by: 0).combine(ch_indel_idx, by: 0).combine(ch_indel_tranches, by: 0)
    GATK4_APPLYVQSR_INDEL ( ch_indel_vqsr_in, fasta, fai, dict, allelespecific, truthsensitivity, 'INDEL' )

    GATK4_SELECTVARIANTS_NORECAL (ch_select_variants_in)

    ch_merge_recal_in = GATK4_SELECTVARIANTS_NORECAL.out.vcf.mix(GATK4_APPLYVQSR_SNP.out.vcf).mix(GATK4_APPLYVQSR_INDEL.out.vcf).groupTuple(by: 0)
    GATK4_MERGEVCFS_RECALIBRATED(ch_merge_recal_in, dict, true)

    emit:
    versions       = ch_versions                                     // channel: [ versions.yml ]
    select_var_snp_vcf     = GATK4_SELECTVARIANTS_SNP.out.vcf
    select_var_snp_tbi     = GATK4_SELECTVARIANTS_SNP.out.tbi
    recal_snp_file     = GATK4_VARIANTRECALIBRATOR_SNP.out.recal // channel: [ val(meta), [ recal ] ]
    recal_snp_index    = GATK4_VARIANTRECALIBRATOR_SNP.out.idx   // channel: [ val(meta), [ idx ] ]
    recal_snp_tranches = GATK4_VARIANTRECALIBRATOR_SNP.out.tranches // channel: [ val(meta), [ tranches ] ]
    vqsr_snp_vcf       = GATK4_APPLYVQSR_SNP.out.vcf             // channel: [ val(meta), [ vcf ] ]
    vqsr_snp_index     = GATK4_APPLYVQSR_SNP.out.tbi             // channel: [ val(meta), [ tbi ] ]

    select_var_indel_vcf     = GATK4_SELECTVARIANTS_INDEL.out.vcf
    select_var_indel_tbi     = GATK4_SELECTVARIANTS_INDEL.out.tbi
    recal_indel_file     = GATK4_VARIANTRECALIBRATOR_INDEL.out.recal // channel: [ val(meta), [ recal ] ]
    recal_indel_index    = GATK4_VARIANTRECALIBRATOR_INDEL.out.idx   // channel: [ val(meta), [ idx ] ]
    recal_indel_tranches = GATK4_VARIANTRECALIBRATOR_INDEL.out.tranches // channel: [ val(meta), [ tranches ] ]
    vqsr_indel_vcf       = GATK4_APPLYVQSR_INDEL.out.vcf             // channel: [ val(meta), [ vcf ] ]
    vqsr_indel_index     = GATK4_APPLYVQSR_INDEL.out.tbi             // channel: [ val(meta), [ tbi ] ]

    select_var_norecal_vcf     = GATK4_SELECTVARIANTS_NORECAL.out.vcf
    select_var_norecal_tbi     = GATK4_SELECTVARIANTS_NORECAL.out.tbi

    merged_recal_vcf     = GATK4_MERGEVCFS_RECALIBRATED.out.vcf
    merged_recal_tbi     = GATK4_MERGEVCFS_RECALIBRATED.out.tbi
}
