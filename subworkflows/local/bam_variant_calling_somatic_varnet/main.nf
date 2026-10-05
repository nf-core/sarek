//
// VarNet tumor-normal variant calling
//
// For all modules here:
// A when clause condition is defined in the conf/modules.config to determine if the module should be run

include { VARNET_FILTER } from '../../../modules/nf-core/varnet/filter'
include { VARNET_PREDICT } from '../../../modules/nf-core/varnet/predict'
include { HTSLIB_BGZIPTABIX as TABIX_VARNET } from '../../../modules/nf-core/htslib/bgziptabix'

workflow BAM_VARIANT_CALLING_SOMATIC_VARNET {
    take:
    bam // channel: [mandatory] [ meta, normal_bam, normal_bai, tumor_bam, tumor_bai ]
    fasta // channel: [mandatory] [ meta, fasta ]
    fasta_fai // channel: [mandatory] [ meta, fasta_fai ]
    intervals // channel: [mandatory] [ meta, bed ] or [ meta, [] ] if no intervals

    main:
    ch_bam = bam.map { meta, normal_bam, normal_bai, tumor_bam, tumor_bai -> [meta, tumor_bam, tumor_bai, normal_bam, normal_bai] }

    VARNET_FILTER(ch_bam, intervals, fasta, fasta_fai)

    VARNET_PREDICT(ch_bam.join(VARNET_FILTER.out.candidates, by: [0]), fasta, fasta_fai)

    // VarNet writes plain gzip, which tabix cannot index, so recompress as BGZF
    TABIX_VARNET(VARNET_PREDICT.out.vcf.map { meta, vcf -> [meta, vcf, [], []] }, 'compress', true, 'vcf')

    vcf_tbi = TABIX_VARNET.out.output.join(TABIX_VARNET.out.index, by: [0])

    // add variantcaller to meta map
    vcf = vcf_tbi.map { meta, vcf_, _tbi -> [meta + [variantcaller: 'varnet'], vcf_] }
    tbi = vcf_tbi.map { meta, _vcf, tbi_ -> [meta + [variantcaller: 'varnet'], tbi_] }

    emit:
    vcf // channel: [ meta, vcf.gz ]
    tbi // channel: [ meta, vcf.gz.tbi ]
}
