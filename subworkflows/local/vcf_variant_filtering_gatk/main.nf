include { GATK4_CNNSCOREVARIANTS      as CNNSCOREVARIANTS               } from '../../../modules/nf-core/gatk4/cnnscorevariants/main'
include { GATK4_FILTERVARIANTTRANCHES as FILTERVARIANTTRANCHES          } from '../../../modules/nf-core/gatk4/filtervarianttranches/main'

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

    ch_versions = Channel.empty()

    //Don't scatter/gather by intervals, because especially for small regions (targeted or WGS), it easily fails with 0 SNPS in region
    cnn_in = vcf.combine(intervals_bed_combined).map{ meta, vcf, tbi, intervals ->
            new_intervals = intervals.simpleName == "no_intervals" ? [] : intervals
            [meta, vcf, tbi, [], new_intervals]
    }

    CNNSCOREVARIANTS(
        cnn_in,
        fasta,
        fasta_fai,
        dict,
        [],
        []
    )

    cnn_out = CNNSCOREVARIANTS.out.vcf.join(CNNSCOREVARIANTS.out.tbi).combine(intervals_bed_combined)
        .map{   meta, cnn_vcf,cnn_tbi, intervals ->
            new_intervals = intervals.simpleName == "no_intervals" ? [] : intervals
            [meta, cnn_vcf, cnn_tbi, new_intervals]
        }

    FILTERVARIANTTRANCHES(
        cnn_out,
        known_sites,
        known_sites_tbi,
        fasta,
        fasta_fai,
        dict
    )

    // Figure out if using intervals or no_intervals
    filtered_vcf = FILTERVARIANTTRANCHES.out.vcf.map{ meta, vcf ->
                                            [[
                                                id:             meta.sample,
                                                num_intervals:  meta.num_intervals,
                                                patient:        meta.patient,
                                                sample:         meta.sample,
                                                sex:            meta.sex,
                                                status:         meta.status,
                                                variantcaller:  "haplotypecaller"
                                            ], vcf]
                                        }

    ch_versions = ch_versions.mix(CNNSCOREVARIANTS.out.versions)
    ch_versions = ch_versions.mix(FILTERVARIANTTRANCHES.out.versions)

    emit:
    versions = ch_versions
    filtered_vcf
}

