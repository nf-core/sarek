//
// VarNet somatic variant calling on tumor-normal pairs and tumor-only samples
//

include { VARNET_FILTER } from '../../../modules/nf-core/varnet/filter'
include { VARNET_PREDICT } from '../../../modules/nf-core/varnet/predict'
include { HTSLIB_BGZIPTABIX as TABIX_VARNET } from '../../../modules/nf-core/htslib/bgziptabix'

workflow BAM_VARIANT_CALLING_SOMATIC_VARNET {
    take:
    bam // channel: [mandatory] [ meta, normal_bam, normal_bai, tumor_bam, tumor_bai ], with [] for the normal in tumor-only mode
    fasta // channel: [mandatory] [ meta, fasta ]
    fasta_fai // channel: [mandatory] [ meta, fasta_fai ]
    intervals_bed_combined // channel: [mandatory] [ bed ] or [] if no intervals
    resources // channel: [mandatory] folder with the tumor-only germline resources, or [] for tumor-normal pairs

    main:
    ch_bam = bam.map { meta, normal_bam, normal_bai, tumor_bam, tumor_bai -> [meta, tumor_bam, tumor_bai, normal_bam, normal_bai] }
    ch_intervals = intervals_bed_combined.map { bed -> [[id: 'intervals'], bed] }

    ch_dbsnp = resources.map { dir -> [[id: 'varnet_dbsnp'], dir ? file("${dir}/GrCh38-dbSNP-build-156-common-variants-0-indexed-CHROMPOS.npz", checkIfExists: true) : []] }
    ch_gnomad = resources.map { dir -> [[id: 'varnet_gnomad'], dir ? file("${dir}/combined_gnomAD_v4.1_genomes_exomes_dict.hg38.chunks.npz", checkIfExists: true) : []] }
    ch_pon = resources.map { dir -> [[id: 'varnet_pon'], dir ? file("${dir}/1000g_pon.hg38.vcf.gz", checkIfExists: true) : []] }

    VARNET_FILTER(ch_bam, ch_intervals, fasta, fasta_fai, ch_dbsnp, ch_gnomad)

    VARNET_PREDICT(ch_bam.join(VARNET_FILTER.out.candidates, by: [0]), fasta, fasta_fai, ch_pon)

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
