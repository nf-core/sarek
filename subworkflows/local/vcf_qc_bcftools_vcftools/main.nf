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

    // `--skip_tools bcftools` to skip bcftools stats, `--skip_tools vcftools` to skip vcftools
    vcf_bcftools = vcf.filter { !(skip_tools && skip_tools.split(',').contains('bcftools')) }
    vcf_vcftools = vcf.filter { !(skip_tools && skip_tools.split(',').contains('vcftools')) }

    BCFTOOLS_STATS(vcf_bcftools.map{ meta, vcf_ -> [ meta, vcf_, [] ] }, [[:],[]], [[:],[]], [[:],[]], [[:],[]], [[:],[]])
    VCFTOOLS_TSTV_COUNT(vcf_vcftools, target_bed, [])
    VCFTOOLS_TSTV_QUAL(vcf_vcftools, target_bed, [])
    VCFTOOLS_SUMMARY(vcf_vcftools, target_bed, [])

    emit:
    bcftools_stats          = BCFTOOLS_STATS.out.stats
    vcftools_tstv_counts    = VCFTOOLS_TSTV_COUNT.out.tstv_count
    vcftools_tstv_qual      = VCFTOOLS_TSTV_QUAL.out.tstv_qual
    vcftools_filter_summary = VCFTOOLS_SUMMARY.out.filter_summary

}
