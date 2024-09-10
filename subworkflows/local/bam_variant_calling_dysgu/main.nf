//
// dysgu variant calling
//
// For all modules here:
// A when clause condition is defined in the conf/modules.config to determine if the module should be run

include { DYSGU } from '../../../modules/nf-core/dysgu/main'

// Seems to be the consensus on upstream modules implementation too
workflow BAM_VARIANT_CALLING_DYSGU {
    take:
    cram          // channel: [mandatory] [ meta, cram, crai ]
    fasta         // channel: [mandatory] [ meta, fasta ]
    fasta_fai     // channel: [mandatory] [ meta, fasta_fai ]

    main:
    versions = Channel.empty()

    DYSGU(cram, fasta, fasta_fai)

    dysgu_vcf = DYSGU.out.vcf

    // Only dysgu SV should get annotated
    // add variantcaller to meta map
    vcf = dysgu_vcf.map { meta, vcf -> [ meta + [ variantcaller:'dysgu' ], vcf ] }

    versions = versions.mix(DYSGU.out.versions)

    emit:
    vcf
    versions
}
