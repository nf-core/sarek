//
// Run Benchmark on called VCF files
//

include { HAPPY_PREPY } from '../../../modules/nf-core/happy/prepy/main'
include { HAPPY_HAPPY } from '../../../modules/nf-core/happy/happy/main'

workflow VCF_BENCHMARK {
    take:
        vcfs
        confidence_bed
        truth_vcf
        fasta
        fasta_fai

    main:

    HAPPY_HAPPY()

    ch_versions = Channel.empty()

    emit:
    versions = ch_versions


}
