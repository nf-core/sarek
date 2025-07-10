//
// TIDDIT single sample variant calling
//
// For all modules here:
// A when clause condition is defined in the conf/modules.config to determine if the module should be run

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

    versions = Channel.empty()

    TIDDIT_NORMAL(cram_normal, fasta, bwa)
    TIDDIT_TUMOR(cram_tumor, fasta, bwa)

    SVDB_MERGE(TIDDIT_NORMAL.out.vcf.join(TIDDIT_TUMOR.out.vcf, failOnDuplicate: true, failOnMismatch: true).map{ meta, vcf_normal, vcf_tumor -> [ meta, [vcf_normal, vcf_tumor] ] }, false, true)

    vcf = SVDB_MERGE.out.vcf

    versions = versions.mix(TIDDIT_NORMAL.out.versions)
    versions = versions.mix(TIDDIT_TUMOR.out.versions)
    versions = versions.mix(SVDB_MERGE.out.versions)

    emit:
    versions
    vcf
}
