include { GATK4_CNNSCOREVARIANTS      } from '../../../modules/nf-core/gatk4/cnnscorevariants/main'
include { GATK4_FILTERVARIANTTRANCHES } from '../../../modules/nf-core/gatk4/filtervarianttranches/main'

workflow VCF_VARIANT_FILTERING_GATK {

    take:
    vcf             // meta, vcf, tbi, intervals
    fasta
    fasta_fai
    dict
    intervals_bed_combined
    known_sites
    known_sites_tbi

    main:

    versions = Channel.empty()

    // Don't scatter/gather by intervals, because especially for small regions (targeted or WGS), it easily fails with 0 SNPS in region
    cnn_in = vcf.combine(intervals_bed_combined).map{ meta, vcf, tbi, intervals -> [ meta, vcf, tbi, [], intervals ] }

    GATK4_CNNSCOREVARIANTS(cnn_in, fasta, fasta_fai, dict, [], [])

    GATK4_FILTERVARIANTTRANCHES(GATK4_CNNSCOREVARIANTS.out.vcf.join(GATK4_CNNSCOREVARIANTS.out.tbi, failOnDuplicate: true, failOnMismatch: true).combine(intervals_bed_combined), known_sites, known_sites_tbi, fasta, fasta_fai, dict)

    filtered_vcf = GATK4_FILTERVARIANTTRANCHES.out.vcf
        // add variantcaller to meta map and remove no longer necessary field: num_intervals
        .map{ meta, vcf -> [ meta - meta.subMap('num_intervals') + [ variantcaller:'haplotypecaller' ], vcf ] }

    versions = versions.mix(GATK4_CNNSCOREVARIANTS.out.versions)
    versions = versions.mix(GATK4_FILTERVARIANTTRANCHES.out.versions)

    emit:
    filtered_vcf

    versions
}
