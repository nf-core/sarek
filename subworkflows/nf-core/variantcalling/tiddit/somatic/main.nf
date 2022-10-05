include { RUN_TIDDIT as RUN_TIDDIT_NORMAL           } from '../single/main.nf'
include { RUN_TIDDIT as RUN_TIDDIT_TUMOR            } from '../single/main.nf'
include { SVDB_MERGE                                } from '../../../../../modules/nf-core/svdb/merge/main.nf'

workflow RUN_TIDDIT_SOMATIC {
    take:
        cram_normal
        cram_tumor
        fasta
        bwa

    main:

    ch_versions = Channel.empty()

    RUN_TIDDIT_NORMAL(cram_normal, fasta, bwa)
    RUN_TIDDIT_TUMOR(cram_tumor, fasta, bwa)

    SVDB_MERGE(RUN_TIDDIT_NORMAL.out.tiddit_vcf.join(RUN_TIDDIT_TUMOR.out.tiddit_vcf)
                                                .map{meta, vcf_normal, vcf_tumor ->
                                                    [meta, [vcf_normal, vcf_tumor]]
                                                }, false)
    tiddit_vcf = SVDB_MERGE.out.vcf

    ch_versions = ch_versions.mix(RUN_TIDDIT_NORMAL.out.versions)
    ch_versions = ch_versions.mix(RUN_TIDDIT_TUMOR.out.versions)
    ch_versions = ch_versions.mix(SVDB_MERGE.out.versions)

    emit:
    versions = ch_versions
    tiddit_vcf
}
