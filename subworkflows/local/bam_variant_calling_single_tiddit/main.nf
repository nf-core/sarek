//
// TIDDIT single sample variant calling
//
// For all modules here:
// A when clause condition is defined in the conf/modules.config to determine if the module should be run

include { HTSLIB_BGZIPTABIX as TABIX_BGZIP_TIDDIT_SV } from '../../../modules/nf-core/htslib/bgziptabix/main'
include { TIDDIT_SV                                 } from '../../../modules/nf-core/tiddit/sv/main'

workflow BAM_VARIANT_CALLING_SINGLE_TIDDIT {
    take:
    cram
    fasta
    fasta_fai
    bwa

    main:

    TIDDIT_SV(cram, fasta.combine(fasta_fai).map { fasta_meta, fasta_path, _fai_meta, fai_path -> [ fasta_meta, fasta_path, fai_path ] }, bwa)

    TABIX_BGZIP_TIDDIT_SV(TIDDIT_SV.out.vcf.map { meta, vcf -> [ meta, vcf, [], [] ] }, 'compress', true, 'vcf')

    ploidy = TIDDIT_SV.out.ploidy
    vcf = TABIX_BGZIP_TIDDIT_SV.out.output.join(TABIX_BGZIP_TIDDIT_SV.out.index).map { meta, gz, _tbi -> [meta + [variantcaller: 'tiddit'], gz] }
    tbi = TABIX_BGZIP_TIDDIT_SV.out.output.join(TABIX_BGZIP_TIDDIT_SV.out.index).map { meta, _gz, tbi -> [meta + [variantcaller: 'tiddit'], tbi] }

    emit:
    ploidy
    vcf
    tbi
}
