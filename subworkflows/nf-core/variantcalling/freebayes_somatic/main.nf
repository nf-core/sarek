include { TABIX_BGZIP as BGZIP_VC_FREEBAYES } from '../../../../modules/nf-core/modules/tabix/bgzip/main'
include { CONCAT_VCF as CONCAT_FREEBAYES    } from '../../../../modules/local/concat_vcf/main'
include { FREEBAYES as FREEBAYES_SOMATIC    } from '../../../../modules/nf-core/modules/freebayes/main'
include { TABIX_TABIX as TABIX_VC_FREEBAYES } from '../../../../modules/nf-core/modules/tabix/tabix/main'

workflow RUN_FREEBAYES_SOMATIC {
    take:
    cram                     // channel: [mandatory] [meta, normalcram, crai, tumor_Cram, tumor_crai, interval]
    fasta                    // channel: [mandatory]
    fasta_fai                // channel: [mandatory]
    intervals_bed_gz         // channel: [optional]  Contains a bed.gz file of all intervals combined provided with the cram input(s). Mandatory if interval files are used.

    main:

    ch_versions = Channel.empty()

    FREEBAYES_SOMATIC(
        cram,
        fasta,
        fasta_fai,
        [], [], [])

    FREEBAYES.out.vcf.branch{
            intervals:    it[0].num_intervals > 1
            no_intervals: it[0].num_intervals <= 1
        }.set{freebayes_vcf_out}

    // Only when no intervals
    TABIX_VC_FREEBAYES(freebayes_vcf_out.no_intervals)

    // Only when using intervals
    BGZIP_VC_FREEBAYES(freebayes_vcf_out.intervals)

    CONCAT_FREEBAYES(
        BGZIP_VC_FREEBAYES.out.output
            .map{ meta, vcf ->
                new_meta = meta.clone()
                new_meta.id = new_meta.sample

                def groupKey = groupKey(meta, meta.num_intervals)
                [new_meta, vcf]
            }.groupTuple(),
        fasta_fai,
        intervals_bed_gz)

    // Mix output channels for "no intervals" and "with intervals" results
    freebayes_vcf = Channel.empty().mix(
                        CONCAT_FREEBAYES.out.vcf,
                        freebayes_vcf_out.no_intervals)
                    .map{ meta, vcf ->
                        meta.variantcaller = "FreeBayes"
                        [meta, vcf]
                    }

    ch_versions = ch_versions.mix(BGZIP_VC_FREEBAYES.out.versions)
    ch_versions = ch_versions.mix(CONCAT_FREEBAYES.out.versions)
    ch_versions = ch_versions.mix(FREEBAYES.out.versions)
    ch_versions = ch_versions.mix(TABIX_VC_FREEBAYES.out.versions)

    emit:
    freebayes_vcf
    versions = ch_versions
}
