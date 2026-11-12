include { BCFTOOLS_STATS                  } from '../../../modules/nf-core/bcftools/stats/main'
include { VCFTOOLS as VCFTOOLS_SUMMARY    } from '../../../modules/nf-core/vcftools/main'
include { VCFTOOLS as VCFTOOLS_TSTV_COUNT } from '../../../modules/nf-core/vcftools/main'
include { VCFTOOLS as VCFTOOLS_TSTV_QUAL  } from '../../../modules/nf-core/vcftools/main'

workflow VCF_QC_BCFTOOLS_VCFTOOLS {
    take:
    vcf
    target_bed
    skip_tools // string: --skip_tools

    main:

    bcftools_stats          = channel.empty()
    vcftools_tstv_counts    = channel.empty()
    vcftools_tstv_qual      = channel.empty()
    vcftools_filter_summary = channel.empty()

    // `--skip_tools bcftools` to skip bcftools stats
    if (!(skip_tools && skip_tools.split(',').contains('bcftools'))) {
        BCFTOOLS_STATS(vcf.map{ meta, vcf_ -> [ meta, vcf_, [] ] }, [[:],[]], [[:],[]], [[:],[]], [[:],[]], [[:],[]])
        bcftools_stats = BCFTOOLS_STATS.out.stats
    }

    // `--skip_tools vcftools` to skip vcftools
    if (!(skip_tools && skip_tools.split(',').contains('vcftools'))) {
        VCFTOOLS_TSTV_COUNT(vcf, target_bed, [])
        VCFTOOLS_TSTV_QUAL(vcf, target_bed, [])
        VCFTOOLS_SUMMARY(vcf, target_bed, [])
        vcftools_tstv_counts    = VCFTOOLS_TSTV_COUNT.out.tstv_count
        vcftools_tstv_qual      = VCFTOOLS_TSTV_QUAL.out.tstv_qual
        vcftools_filter_summary = VCFTOOLS_SUMMARY.out.filter_summary
    }

    emit:
    bcftools_stats
    vcftools_tstv_counts
    vcftools_tstv_qual
    vcftools_filter_summary

}
