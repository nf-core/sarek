//
// JOINT GERMLINE CALLING OR SINGLE-SAMPLE FILTERING
//
// For callers producing GATK-compatible gVCFs and VCFs (GATK and Parabricks HaplotypeCaller):
// with joint_germline, joint-genotype the gVCFs; otherwise CNN-filter the single-sample VCFs, unless skip_filter
//

include { BAM_JOINT_CALLING_GERMLINE_GATK       } from '../bam_joint_calling_germline_gatk/main'
include { VCF_VARIANT_FILTERING_HAPLOTYPECALLER } from '../vcf_variant_filtering_haplotypecaller/main'

workflow VCF_JOINT_CALLING_OR_FILTERING_HAPLOTYPECALLER {
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
    variant_caller         // string:  name of the caller, used in meta and output paths of the joint calls
    joint_germline         // boolean: joint-genotype the gVCFs instead of filtering the single-sample VCFs
    skip_filter            // boolean: pass the single-sample VCFs on unfiltered

    main:
    BAM_JOINT_CALLING_GERMLINE_GATK(
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

    VCF_VARIANT_FILTERING_HAPLOTYPECALLER(
        vcf.join(tbi, failOnDuplicate: true, failOnMismatch: true).filter { _meta, _vcf, _tbi -> !joint_germline && !skip_filter },
        fasta.map{ _meta, fasta_ -> [ fasta_ ] },
        fasta_fai.map{ _meta, fasta_fai_ -> [ fasta_fai_ ] },
        dict.map{ _meta, dict_ -> [ dict_ ] },
        intervals_bed_combined,
        known_sites_indels.concat(known_sites_snps).flatten().unique().collect(),
        known_sites_indels_tbi.concat(known_sites_snps_tbi).flatten().unique().collect())

    emit:
    vcf = BAM_JOINT_CALLING_GERMLINE_GATK.out.genotype_vcf
        .mix(VCF_VARIANT_FILTERING_HAPLOTYPECALLER.out.filtered_vcf)
        .mix(vcf.filter { _meta, _vcf -> !joint_germline && skip_filter })
    tbi = BAM_JOINT_CALLING_GERMLINE_GATK.out.genotype_index
        .mix(VCF_VARIANT_FILTERING_HAPLOTYPECALLER.out.filtered_tbi)
        .mix(tbi.filter { _meta, _tbi -> !joint_germline && skip_filter })
}
