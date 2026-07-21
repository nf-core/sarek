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
    versions = channel.empty()
    generate_pon = false

    // cram carries [ meta, tumor, normal ]; the module now expects tumor/normal index
    // slots too. The incoming channel does not carry indexes, so pass [] placeholders
    // (the script converts/indexes CRAMs internally and never reads these inputs).
    cram_input = cram.map { meta, tumor, normal -> [meta, tumor, [], normal, []] }

    // module now takes fasta and fasta_fai as a single [ meta, fasta, fasta_fai ] tuple
    fasta_input = fasta.combine(fasta_fai).map { meta, fasta_, _meta_fai, fasta_fai_ -> [meta, fasta_, fasta_fai_] }

    CNVKIT_BATCH(cram_input, fasta_input, targets, reference, generate_pon)

    // right now we do not use an input VCF to improve the calling of B alleles
    // based on SNV frequencies from the VCF file
    // in the future we might consider to add this, by connecting the emission from
    // SNV variant calling modules
    CNVKIT_CALL(CNVKIT_BATCH.out.cns.map{ meta, cns -> [meta, cns[2], []]})

    // export to VCF for compatibility with other tools
    CNVKIT_EXPORT(CNVKIT_CALL.out.cns)

    ch_genemetrics = CNVKIT_BATCH.out.cnr.join(CNVKIT_BATCH.out.cns).map{ meta, cnr, cns -> [meta, cnr, cns[2]]}
    CNVKIT_GENEMETRICS(ch_genemetrics)

    emit:
    cnv_calls_raw    = CNVKIT_CALL.out.cns      // channel: [ meta, cns ]
    cnv_calls_export = CNVKIT_EXPORT.out.output // channel: [ meta, export_format ]
    versions                                    // channel: [ versions.yml ]
}
