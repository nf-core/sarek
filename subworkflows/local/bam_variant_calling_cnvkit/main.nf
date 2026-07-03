//
// CNVKIT calling
//
// For all modules here:
// A when clause condition is defined in the conf/modules.config to determine if the module should be run

include { CNVKIT_BATCH       } from '../../../modules/nf-core/cnvkit/batch/main'
include { CNVKIT_CALL        } from '../../../modules/nf-core/cnvkit/call/main'
include { CNVKIT_EXPORT      } from '../../../modules/nf-core/cnvkit/export/main'
include { CNVKIT_GENEMETRICS } from '../../../modules/nf-core/cnvkit/genemetrics/main'

workflow BAM_VARIANT_CALLING_CNVKIT {
    take:
    cram                // channel: [mandatory] meta, cram
    fasta               // channel: [mandatory] meta, fasta
    fasta_fai           // channel: [optional]  meta, fasta_fai
    targets             // channel: [mandatory] meta, bed
    reference           // channel: [optional]  meta, cnn

    main:
    versions = Channel.empty()
    generate_pon = false

    CNVKIT_BATCH(cram, fasta, fasta_fai, targets, reference, generate_pon)

    // right now we do not use an input VCF to improve the calling of B alleles
    // based on SNV frequencies from the VCF file
    // in the future we might consider to add this, by connecting the emission from
    // SNV variant calling modules
    CNVKIT_CALL(CNVKIT_BATCH.out.cns.map{ meta, cns -> [meta, cns[2], []]})

    // export to VCF for compatibility with other tools
    CNVKIT_EXPORT(CNVKIT_CALL.out.cns)

    ch_genemetrics = CNVKIT_BATCH.out.cnr.join(CNVKIT_BATCH.out.cns).map{ meta, cnr, cns -> [meta, cnr, cns[2]]}
    CNVKIT_GENEMETRICS(ch_genemetrics)

    versions = versions.mix(CNVKIT_BATCH.out.versions)
    versions = versions.mix(CNVKIT_GENEMETRICS.out.versions)
    versions = versions.mix(CNVKIT_CALL.out.versions)
    versions = versions.mix(CNVKIT_EXPORT.out.versions)
    emit:
    cnv_calls_raw    = CNVKIT_CALL.out.cns      // channel: [ meta, cns ]
    cnv_calls_export = CNVKIT_EXPORT.out.output // channel: [ meta, export_format ]
    versions                                    // channel: [ versions.yml ]
}
