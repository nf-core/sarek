//
// VarNet somatic variant calling on tumor-normal pairs and tumor-only samples
//

include { VARNET_FILTER } from '../../../modules/nf-core/varnet/filter'
include { VARNET_PREDICT } from '../../../modules/nf-core/varnet/predict'
include { HTSLIB_BGZIPTABIX as TABIX_VARNET } from '../../../modules/nf-core/htslib/bgziptabix'

workflow BAM_VARIANT_CALLING_VARNET {
    take:
    bam // channel: [mandatory] [ meta, normal_bam, normal_bai, tumor_bam, tumor_bai ], with [] for the normal in tumor-only mode
    fasta // channel: [mandatory] [ meta, fasta ]
    fasta_fai // channel: [mandatory] [ meta, fasta_fai ]
    intervals_bed_combined // channel: [mandatory] [ bed ] or [] if no intervals
    dbsnp // channel: [mandatory] [ meta, dbsnp.npz ], with [] for tumor-normal pairs
    gnomad // channel: [mandatory] [ meta, gnomad.npz ], with [] for tumor-normal pairs
    pon // channel: [mandatory] [ meta, pon.vcf.gz ], with [] for tumor-normal pairs

    main:
    ch_bam = bam.map { meta, normal_bam, normal_bai, tumor_bam, tumor_bai -> [meta, tumor_bam, tumor_bai, normal_bam, normal_bai] }
    ch_intervals = intervals_bed_combined.map { bed -> [[id: 'intervals'], bed] }

    VARNET_FILTER(ch_bam, ch_intervals, fasta, fasta_fai, dbsnp, gnomad)

    VARNET_PREDICT(ch_bam.join(VARNET_FILTER.out.candidates, by: [0]), fasta, fasta_fai, pon)

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
