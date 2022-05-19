include { TABIX_BGZIP as BGZIP_VC_DEEPVARIANT_GVCF } from '../../../../modules/nf-core/modules/tabix/bgzip/main'
include { TABIX_BGZIP as BGZIP_VC_DEEPVARIANT_VCF  } from '../../../../modules/nf-core/modules/tabix/bgzip/main'
include { CONCAT_VCF as CONCAT_DEEPVARIANT_GVCF    } from '../../../../modules/local/concat_vcf/main'
include { CONCAT_VCF as CONCAT_DEEPVARIANT_VCF     } from '../../../../modules/local/concat_vcf/main'
include { DEEPVARIANT                              } from '../../../../modules/nf-core/modules/deepvariant/main'
include { TABIX_TABIX as TABIX_VC_DEEPVARIANT_GVCF } from '../../../../modules/nf-core/modules/tabix/tabix/main'
include { TABIX_TABIX as TABIX_VC_DEEPVARIANT_VCF  } from '../../../../modules/nf-core/modules/tabix/tabix/main'

//TODO: benchmark if it is better to provide multiple bed files & run on multiple machines + mergeing afterwards || one containing all intervals and run on one larger machine
// Deepvariant: https://github.com/google/deepvariant/issues/510
workflow RUN_DEEPVARIANT {
    take:
    cram                     // channel: [mandatory] [meta, cram, crai, interval]
    fasta                    // channel: [mandatory]
    fasta_fai                // channel: [mandatory]
    intervals_bed_gz         // channel: [optional]  Contains a bed.gz file of all intervals combined provided with the cram input(s). Mandatory if interval files are used.

    main:

    ch_versions = Channel.empty()

    DEEPVARIANT(cram, fasta, fasta_fai)

    DEEPVARIANT.out.vcf.branch{
        intervals:    it[0].num_intervals > 1
        no_intervals: it[0].num_intervals <= 1
    }.set{deepvariant_vcf_out}

    DEEPVARIANT.out.gvcf.branch{
        intervals:    it[0].num_intervals > 1
        no_intervals: it[0].num_intervals <= 1
    }.set{deepvariant_gvcf_out}

    // Only when no intervals
    TABIX_VC_DEEPVARIANT_VCF(deepvariant_vcf_out.no_intervals)
    TABIX_VC_DEEPVARIANT_GVCF(deepvariant_gvcf_out.no_intervals)

    // Only when using intervals
    BGZIP_VC_DEEPVARIANT_VCF(deepvariant_vcf_out.intervals)
    BGZIP_VC_DEEPVARIANT_GVCF(deepvariant_gvcf_out.intervals)

    CONCAT_DEEPVARIANT_VCF(
        BGZIP_VC_DEEPVARIANT_VCF.out.output
            .map{ meta, vcf ->

                new_meta = [patient:meta.patient, sample:meta.sample, status:meta.status, gender:meta.gender, id:meta.sample, num_intervals:meta.num_intervals]

                def groupKey = groupKey(new_meta, meta.num_intervals)
                [new_meta, vcf]
            }.groupTuple(),
        fasta_fai,
        intervals_bed_gz)

    CONCAT_DEEPVARIANT_GVCF(
        BGZIP_VC_DEEPVARIANT_GVCF.out.output
            .map{ meta, vcf ->

                new_meta = [patient:meta.patient, sample:meta.sample, status:meta.status, gender:meta.gender, id:meta.sample, num_intervals:meta.num_intervals]

                def groupKey = groupKey(new_meta, meta.num_intervals)
                [new_meta, vcf]
            }.groupTuple(),
        fasta_fai,
        intervals_bed_gz)

    // Mix output channels for "no intervals" and "with intervals" results
    deepvariant_vcf = Channel.empty().mix(
                        CONCAT_DEEPVARIANT_GVCF.out.vcf,
                        CONCAT_DEEPVARIANT_VCF.out.vcf,
                        deepvariant_gvcf_out.no_intervals,
                        deepvariant_vcf_out.no_intervals)
                    .map{ meta, vcf ->
                        meta.variantcaller = "Deepvariant"
                        [meta, vcf]
                    }

    ch_versions = ch_versions.mix(BGZIP_VC_DEEPVARIANT_GVCF.out.versions)
    ch_versions = ch_versions.mix(BGZIP_VC_DEEPVARIANT_VCF.out.versions)
    ch_versions = ch_versions.mix(CONCAT_DEEPVARIANT_GVCF.out.versions)
    ch_versions = ch_versions.mix(CONCAT_DEEPVARIANT_VCF.out.versions)
    ch_versions = ch_versions.mix(DEEPVARIANT.out.versions)
    ch_versions = ch_versions.mix(TABIX_VC_DEEPVARIANT_GVCF.out.versions)
    ch_versions = ch_versions.mix(TABIX_VC_DEEPVARIANT_VCF.out.versions)

    emit:
    deepvariant_vcf
    versions = ch_versions
}
