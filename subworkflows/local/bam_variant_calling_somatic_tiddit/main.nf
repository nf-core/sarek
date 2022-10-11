include { BAM_VARIANT_CALLING_SINGLE_TIDDIT as TIDDIT_NORMAL } from '../bam_variant_calling_single_tiddit/main.nf'
include { BAM_VARIANT_CALLING_SINGLE_TIDDIT as TIDDIT_TUMOR  } from '../bam_variant_calling_single_tiddit/main.nf'
include { SVDB_MERGE                                         } from '../../../modules/nf-core/svdb/merge/main.nf'

workflow BAM_VARIANT_CALLING_SOMATIC_TIDDIT {
    take:
        cram_normal
        cram_tumor
        fasta
        bwa

    main:

    ch_versions = Channel.empty()

    TIDDIT_NORMAL(cram_normal, fasta, bwa)
    TIDDIT_TUMOR(cram_tumor, fasta, bwa)

    SVDB_MERGE(TIDDIT_NORMAL.out.tiddit_vcf.join(TIDDIT_TUMOR.out.tiddit_vcf)
                                                .map{meta, vcf_normal, vcf_tumor ->
                                                    [meta, [vcf_normal, vcf_tumor]]
                                                }, false)
    tiddit_vcf = SVDB_MERGE.out.vcf

    ch_versions = ch_versions.mix(TIDDIT_NORMAL.out.versions)
    ch_versions = ch_versions.mix(TIDDIT_TUMOR.out.versions)
    ch_versions = ch_versions.mix(SVDB_MERGE.out.versions)

    emit:
    versions = ch_versions
    tiddit_vcf
}
