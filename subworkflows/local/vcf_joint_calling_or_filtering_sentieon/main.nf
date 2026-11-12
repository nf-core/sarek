//
// JOINT GERMLINE CALLING OR SINGLE-SAMPLE FILTERING (SENTIEON)
//
// For Sentieon DnaScope and Haplotyper: with joint_germline, joint-genotype the gVCFs with Sentieon GVCFtyper;
// otherwise filter the single-sample VCFs (DnaScope: Sentieon DNAModelApply, Haplotyper: CNN filtering), unless skip_filter
//

include { BAM_JOINT_CALLING_GERMLINE_SENTIEON   } from '../bam_joint_calling_germline_sentieon/main'
include { SENTIEON_DNAMODELAPPLY                } from '../../../modules/nf-core/sentieon/dnamodelapply/main'
include { VCF_VARIANT_FILTERING_HAPLOTYPECALLER } from '../vcf_variant_filtering_haplotypecaller/main'

workflow VCF_JOINT_CALLING_OR_FILTERING_SENTIEON {
    take:
    gvcf_tbi_intervals     // channel: [ meta, gvcf, tbi, intervals ]
    vcf                    // channel: [ meta, vcf ]
    tbi                    // channel: [ meta, tbi ]
    fasta                  // channel: [ meta, fasta ]
    fasta_fai              // channel: [ meta, fasta_fai ]
    dict                   // channel: [ meta, dict ]
    dbsnp
    dbsnp_tbi
    dbsnp_vqsr
    known_sites_indels
    known_sites_indels_tbi
    known_indels_vqsr
    known_sites_snps
    known_sites_snps_tbi
    known_snps_vqsr
    intervals_bed_combined // channel: intervals/target regions in one file unzipped, no_intervals.bed if no_intervals
    dnascope_model         // channel: [ model ], only used for sentieon_dnascope
    variant_caller         // string:  'sentieon_dnascope' or 'sentieon_haplotyper'
    joint_germline         // boolean: joint-genotype the gVCFs instead of filtering the single-sample VCFs
    skip_filter            // boolean: pass the single-sample VCFs on unfiltered

    main:
    BAM_JOINT_CALLING_GERMLINE_SENTIEON(
        gvcf_tbi_intervals.filter { _meta, _gvcf, _tbi, _intervals -> joint_germline },
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
        known_snps_vqsr,
        variant_caller)

    vcf_tbi_to_filter = vcf.join(tbi, failOnDuplicate: true, failOnMismatch: true).filter { _meta, _vcf, _tbi -> !joint_germline && !skip_filter }

    SENTIEON_DNAMODELAPPLY(
        vcf_tbi_to_filter.filter { _meta, _vcf, _tbi -> variant_caller == 'sentieon_dnascope' },
        fasta,
        fasta_fai,
        dnascope_model.map{ model -> [ [ id:model.baseName ], model ] })

    VCF_VARIANT_FILTERING_HAPLOTYPECALLER(
        vcf_tbi_to_filter.filter { _meta, _vcf, _tbi -> variant_caller == 'sentieon_haplotyper' },
        fasta.map{ _meta, fasta_ -> [ fasta_ ] },
        fasta_fai.map{ _meta, fasta_fai_ -> [ fasta_fai_ ] },
        dict.map{ _meta, dict_ -> [ dict_ ] },
        intervals_bed_combined,
        known_sites_indels.concat(known_sites_snps).flatten().unique().collect(),
        known_sites_indels_tbi.concat(known_sites_snps_tbi).flatten().unique().collect())

    emit:
    vcf = BAM_JOINT_CALLING_GERMLINE_SENTIEON.out.genotype_vcf
        .mix(SENTIEON_DNAMODELAPPLY.out.vcf)
        .mix(VCF_VARIANT_FILTERING_HAPLOTYPECALLER.out.filtered_vcf)
        .mix(vcf.filter { _meta, _vcf -> !joint_germline && skip_filter })
    tbi = BAM_JOINT_CALLING_GERMLINE_SENTIEON.out.genotype_index
        .mix(SENTIEON_DNAMODELAPPLY.out.tbi)
        .mix(VCF_VARIANT_FILTERING_HAPLOTYPECALLER.out.filtered_tbi)
        .mix(tbi.filter { _meta, _tbi -> !joint_germline && skip_filter })
}
