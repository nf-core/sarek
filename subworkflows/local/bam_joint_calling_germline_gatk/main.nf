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
include { TABIX_TABIX as TABIX                                   } from '../../../modules/nf-core/tabix/tabix/main'

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
        .map{ meta, gvcf, tbi, intervals -> [ [ intervals_name:intervals.simpleName, id:'joint_variant_calling', num_intervals:meta.num_intervals ], gvcf, tbi, intervals ]}
        .groupTuple(by:[0, 3])
        .map{ meta, gvcf, tbi, intervals -> [ meta, gvcf, tbi, intervals, [], [] ] }

    // Convert all sample vcfs into a genomicsdb workspace using genomicsdbimport
    GATK4_GENOMICSDBIMPORT ( gendb_input, false, false, false )

    genotype_input = GATK4_GENOMICSDBIMPORT.out.genomicsdb.map{ meta, genomicsdb -> [ meta, genomicsdb, [], [], [] ] }

    // Joint genotyping performed using GenotypeGVCFs
    // Sort vcfs called by interval within each VCF

    GATK4_GENOTYPEGVCFS (genotype_input, fasta, fai, dict, dbsnp, dbsnp_tbi)

    BCFTOOLS_SORT(GATK4_GENOTYPEGVCFS.out.vcf)
    vcfs_sorted_input = BCFTOOLS_SORT.out.vcf.branch{
        intervals:    it[0].num_intervals > 1
        no_intervals: it[0].num_intervals <= 1
    }

    vcfs_sorted_input_no_intervals =  vcfs_sorted_input.no_intervals.map{ meta, vcf ->
        [ meta.subMap('num_intervals') + [ id:'joint_variant_calling', patient:'all_samples', variantcaller:'haplotypecaller' ], vcf ]
    }

    // Index vcf files if no scatter/gather by intervals
    TABIX(vcfs_sorted_input_no_intervals)

    // Merge scatter/gather vcfs & index
    // Rework meta for variantscalled.csv and annotation tools
    MERGE_GENOTYPEGVCFS(
        vcfs_sorted_input.intervals.map{ meta, vcf ->
            [ meta.subMap('num_intervals') + [ id:'joint_variant_calling', patient:'all_samples', variantcaller:'haplotypecaller' ], vcf ]
        }.groupTuple(),
        dict.map{ it -> [ [ id:'dict' ], it ] } )

    vqsr_input = Channel.empty().mix(
        MERGE_GENOTYPEGVCFS.out.vcf.join(MERGE_GENOTYPEGVCFS.out.tbi),
        vcfs_sorted_input_no_intervals.join(TABIX.out.tbi)
    )

    // Group resource labels for SNP and INDEL
    snp_resource_labels   = Channel.empty().mix(known_snps_vqsr,dbsnp_vqsr).collect()
    indel_resource_labels = Channel.empty().mix(known_indels_vqsr,dbsnp_vqsr).collect()

    // Recalibrate SNP and INDEL separately.
    VARIANTRECALIBRATOR_SNP(
        vqsr_input,
        resource_snps_vcf,
        resource_snps_tbi,
        snp_resource_labels,
        fasta,
        fai,
        dict)

    VARIANTRECALIBRATOR_INDEL(
        vqsr_input,
        resource_indels_vcf,
        resource_indels_tbi,
        indel_resource_labels,
        fasta,
        fai,
        dict)

    //Prepare SNP and INDEL separately for ApplyVQSR

    // Join results of variant recalibration into a single channel tuple
    // Rework meta for variantscalled.csv and annotation tools
    vqsr_input_snp   = vqsr_input.join(VARIANTRECALIBRATOR_SNP.out.recal)
        .join(VARIANTRECALIBRATOR_SNP.out.idx)
        .join(VARIANTRECALIBRATOR_SNP.out.tranches)
        .map{ meta, vcf, tbi, recal, index, tranche ->
            [ meta.subMap('num_intervals') + [ id:'recalibrated_joint_variant_calling', patient:'all_samples', variantcaller:'haplotypecaller' ], vcf, tbi, recal, index, tranche ]
        }

    // Join results of variant recalibration into a single channel tuple
    // Rework meta for variantscalled.csv and annotation tools
    vqsr_input_indel = vqsr_input.join(VARIANTRECALIBRATOR_INDEL.out.recal)
        .join(VARIANTRECALIBRATOR_INDEL.out.idx)
        .join(VARIANTRECALIBRATOR_INDEL.out.tranches)
        .map{ meta, vcf, tbi, recal, index, tranche ->
            [ meta.subMap('num_intervals') + [ id:'recalibrated_joint_variant_calling', patient:'all_samples', variantcaller:'haplotypecaller' ], vcf, tbi, recal, index, tranche ]
        }

    GATK4_APPLYVQSR_SNP(
        vqsr_input_snp,
        fasta,
        fai,
        dict)

    GATK4_APPLYVQSR_INDEL(
        vqsr_input_indel,
        fasta,
        fai,
        dict)

    vqsr_snp_vcf = GATK4_APPLYVQSR_SNP.out.vcf
    vqsr_indel_vcf = GATK4_APPLYVQSR_INDEL.out.vcf

    //Merge VQSR outputs into final VCF
    MERGE_VQSR(
        vqsr_snp_vcf.mix(vqsr_indel_vcf).groupTuple(),
        dict.map{ it -> [ [ id:'dict' ], it ] }
    )

    genotype_vcf   = Channel.empty().mix(vcfs_sorted_input_no_intervals, MERGE_GENOTYPEGVCFS.out.vcf, MERGE_VQSR.out.vcf)
    genotype_index = Channel.empty().mix(TABIX.out.tbi, MERGE_GENOTYPEGVCFS.out.tbi, MERGE_VQSR.out.tbi)

    versions = versions.mix(GATK4_GENOMICSDBIMPORT.out.versions)
    versions = versions.mix(GATK4_GENOTYPEGVCFS.out.versions)
    versions = versions.mix(VARIANTRECALIBRATOR_SNP.out.versions)
    versions = versions.mix(GATK4_APPLYVQSR_SNP.out.versions)

    emit:
    genotype_index  // channel: [ val(meta), [ tbi ] ]
    genotype_vcf    // channel: [ val(meta), [ vcf ] ]

    versions        // channel: [ versions.yml ]
}
