//
// JOINT GERMLINE CALLING
//
// Merge samples with genomicsdbimport, perform joint genotyping with genotypeGVCFS
//

include { BCFTOOLS_SORT                                          } from '../../../modules/nf-core/bcftools/sort/main'
include { GATK4_APPLYVQSR           as GATK4_APPLYVQSR_INDEL     } from '../../../modules/nf-core/gatk4/applyvqsr/main'
include { GATK4_APPLYVQSR           as GATK4_APPLYVQSR_SNP       } from '../../../modules/nf-core/gatk4/applyvqsr/main'
include { GATK4_GENOMICSDBIMPORT                                 } from '../../../modules/nf-core/gatk4/genomicsdbimport/main'
include { GATK4_GENOTYPEGVCFS                                    } from '../../../modules/nf-core/gatk4/genotypegvcfs/main'
include { GATK4_MERGEVCFS           as MERGE_GENOTYPEGVCFS       } from '../../../modules/nf-core/gatk4/mergevcfs/main'
include { GATK4_MERGEVCFS           as MERGE_VQSR                } from '../../../modules/nf-core/gatk4/mergevcfs/main'
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

    main:
    versions = Channel.empty()

    // Map input for GenomicsDBImport
    // Rename based on num_intervals, group all samples by their interval_name/interval_file and restructure for channel
    // Group by [0, 3] to avoid a list of metas and make sure that any intervals
    gendb_input = input
        .map{ meta, gvcf, tbi, intervals -> [ [ id:'joint_variant_calling', intervals_name:intervals.baseName, num_intervals:meta.num_intervals ], gvcf, tbi, intervals ] }
        .groupTuple(by:3) //join on interval file
        .map{ meta_list, gvcf, tbi, intervals ->
            // meta is now a list of [meta1, meta2] but they are all the same. So take the first element.
            [ meta_list[0], gvcf, tbi, intervals, [], [] ]
        }

    // Convert all sample vcfs into a genomicsdb workspace using genomicsdbimport
    GATK4_GENOMICSDBIMPORT(gendb_input, false, false, false)

    genotype_input = GATK4_GENOMICSDBIMPORT.out.genomicsdb.map{ meta, genomicsdb -> [ meta, genomicsdb, [], [], [] ] }

    // Joint genotyping performed using GenotypeGVCFs
    // Sort vcfs called by interval within each VCF

    GATK4_GENOTYPEGVCFS(genotype_input, fasta, fai, dict, dbsnp.map{ it -> [ [:], it ] }, dbsnp_tbi.map{ it -> [ [:], it ] })

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
        fasta.map{ meta, fasta -> [ fasta ] },
        fai.map{ meta, fai -> [ fai ] },
        dict.map{ meta, dict -> [ dict ] })

    VARIANTRECALIBRATOR_SNP(
        vqsr_input,
        resource_snps_vcf,
        resource_snps_tbi,
        snps_resource_label,
        fasta.map{ meta, fasta -> [ fasta ] },
        fai.map{ meta, fai -> [ fai ] },
        dict.map{ meta, dict -> [ dict ] })

    //Prepare SNPs and INDELs for ApplyVQSR
    // Step 1. : ApplyVQSR to SNPs
    // Step 2. : Use ApplyVQSR_SNP output and run ApplyVQSR_INDEL. This avoids duplicate entries in the vcf as described here: https://hpc.nih.gov/training/gatk_tutorial/vqsr.html

    // Join results of variant recalibration into a single channel tuple
    // Rework meta for variantscalled.csv and annotation tools
    vqsr_input_snp = vqsr_input.join(VARIANTRECALIBRATOR_SNP.out.recal, failOnDuplicate: true)
        .join(VARIANTRECALIBRATOR_SNP.out.idx, failOnDuplicate: true)
        .join(VARIANTRECALIBRATOR_SNP.out.tranches, failOnDuplicate: true)
        .map{ meta, vcf, tbi, recal, index, tranche -> [ meta + [ id:'recalibrated_joint_variant_calling' ], vcf, tbi, recal, index, tranche ] }

    GATK4_APPLYVQSR_SNP(
        vqsr_input_snp,
        fasta.map{ meta, fasta -> [ fasta ] },
        fai.map{ meta, fai -> [ fai ] },
        dict.map{ meta, dict -> [ dict ] })

    // Join results of ApplyVQSR_SNP and use as input for Indels to avoid duplicate entries in the result
    // Rework meta for variantscalled.csv and annotation tools
    vqsr_input_indel = GATK4_APPLYVQSR_SNP.out.vcf.join(GATK4_APPLYVQSR_SNP.out.tbi).map{ meta, vcf, tbi -> [ meta + [ id:'joint_variant_calling' ], vcf, tbi ]}
        .join(VARIANTRECALIBRATOR_INDEL.out.recal, failOnDuplicate: true)
        .join(VARIANTRECALIBRATOR_INDEL.out.idx, failOnDuplicate: true)
        .join(VARIANTRECALIBRATOR_INDEL.out.tranches, failOnDuplicate: true)
        .map{ meta, vcf, tbi, recal, index, tranche -> [ meta + [ id:'recalibrated_joint_variant_calling' ], vcf, tbi, recal, index, tranche ] }

    GATK4_APPLYVQSR_INDEL(
        vqsr_input_indel,
        fasta.map{ meta, fasta -> [ fasta ] },
        fai.map{ meta, fai -> [ fai ] },
        dict.map{ meta, dict -> [ dict ] })


    // The following is an ugly monster to achieve the following:
    // When MERGE_GENOTYPEGVCFS and GATK4_APPLYVQSR are run, then use output from APPLYVQSR
    // When MERGE_GENOTYPEGVCFS and NOT GATK4_APPLYVQSR , then use the output from MERGE_GENOTYPEGVCFS

    merge_vcf_for_join = MERGE_GENOTYPEGVCFS.out.vcf.map{meta, vcf -> [[id: 'joint_variant_calling'] , vcf]}
    merge_tbi_for_join = MERGE_GENOTYPEGVCFS.out.tbi.map{meta, tbi -> [[id: 'joint_variant_calling'] , tbi]}

    // Remap for both to have the same key, if ApplyBQSR is not run, the channel is empty --> populate with empty elements
    vqsr_vcf_for_join = GATK4_APPLYVQSR_INDEL.out.vcf.ifEmpty([[:], []]).map{meta, vcf -> [[id: 'joint_variant_calling'] , vcf]}
    vqsr_tbi_for_join = GATK4_APPLYVQSR_INDEL.out.tbi.ifEmpty([[:], []]).map{meta, tbi -> [[id: 'joint_variant_calling'] , tbi]}

    // Join on metamap
    // If both --> meta, vcf_merged, vcf_bqsr
    // If not VQSR --> meta, vcf_merged, []
    // if the second is empty, use the first
    genotype_vcf = merge_vcf_for_join.join(vqsr_vcf_for_join, remainder: true).map{
        meta, joint_vcf, recal_vcf ->

        vcf_out = recal_vcf ?: joint_vcf

        [[id:"joint_variant_calling", patient:"all_samples", variantcaller:"haplotypecaller"], vcf_out]
    }

    genotype_index = merge_tbi_for_join.join(vqsr_tbi_for_join, remainder: true).map{
        meta, joint_tbi, recal_tbi ->

        tbi_out = recal_tbi ?: joint_tbi

        [[id:"joint_variant_calling", patient:"all_samples", variantcaller:"haplotypecaller"], tbi_out]
    }

    versions = versions.mix(GATK4_GENOMICSDBIMPORT.out.versions)
    versions = versions.mix(GATK4_GENOTYPEGVCFS.out.versions)
    versions = versions.mix(VARIANTRECALIBRATOR_SNP.out.versions)
    versions = versions.mix(GATK4_APPLYVQSR_SNP.out.versions)

    emit:
    genotype_index  // channel: [ val(meta), [ tbi ] ]
    genotype_vcf    // channel: [ val(meta), [ vcf ] ]

    versions        // channel: [ versions.yml ]
}
