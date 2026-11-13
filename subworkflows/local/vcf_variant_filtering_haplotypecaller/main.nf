include { GATK4_CNNSCOREVARIANTS      as CNNSCOREVARIANTS               } from '../../../modules/nf-core/gatk4/cnnscorevariants/main'
include { GATK4_FILTERVARIANTTRANCHES as FILTERVARIANTTRANCHES          } from '../../../modules/nf-core/gatk4/filtervarianttranches/main'

workflow VCF_VARIANT_FILTERING_HAPLOTYPECALLER {

    take:
    vcf                    // channel: [ meta, vcf, tbi ]
    fasta                  // channel: [ meta, fasta ]
    fasta_fai              // channel: [ meta, fasta_fai ]
    dict                   // channel: [ meta, dict ]
    intervals_bed_combined // channel: intervals/target regions in one file unzipped, no_intervals.bed if no_intervals
    known_sites_indels
    known_sites_indels_tbi
    known_sites_snps
    known_sites_snps_tbi

    main:
    fasta_     = fasta.map{ _meta, fasta_file -> [ fasta_file ] }
    fasta_fai_ = fasta_fai.map{ _meta, fasta_fai_file -> [ fasta_fai_file ] }
    dict_      = dict.map{ _meta, dict_file -> [ dict_file ] }
    // FilterVariantTranches uses every known site as a resource
    known_sites     = known_sites_indels.concat(known_sites_snps).flatten().unique().collect()
    known_sites_tbi = known_sites_indels_tbi.concat(known_sites_snps_tbi).flatten().unique().collect()

    // Don't scatter/gather by intervals, because especially for small regions (targeted or WGS), it easily fails with 0 SNPS in region
    cnn_in = vcf.combine(intervals_bed_combined).map{ meta, vcf_, tbi, intervals -> [ meta, vcf_, tbi, [], intervals ] }

    CNNSCOREVARIANTS(cnn_in, fasta_, fasta_fai_, dict_, [], [])

    FILTERVARIANTTRANCHES(CNNSCOREVARIANTS.out.vcf.join(CNNSCOREVARIANTS.out.tbi, failOnDuplicate: true, failOnMismatch: true).combine(intervals_bed_combined), known_sites, known_sites_tbi, fasta_, fasta_fai_, dict_)

    filtered_vcf = FILTERVARIANTTRANCHES.out.vcf
        // remove no longer necessary field: num_intervals
        .map{ meta, vcf_ -> [ meta - meta.subMap('num_intervals'), vcf_ ] }

    filtered_tbi = FILTERVARIANTTRANCHES.out.tbi
        // remove no longer necessary field: num_intervals
        .map{ meta, tbi -> [ meta - meta.subMap('num_intervals'), tbi ] }


    emit:
    filtered_vcf
    filtered_tbi

}
