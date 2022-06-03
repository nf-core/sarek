include { GATK4_MERGEVCFS as MERGE_DEEPVARIANT_GVCF } from '../../../../modules/nf-core/modules/gatk4/mergevcfs/main'
include { GATK4_MERGEVCFS as MERGE_DEEPVARIANT_VCF  } from '../../../../modules/nf-core/modules/gatk4/mergevcfs/main'
include { DEEPVARIANT                               } from '../../../../modules/nf-core/modules/deepvariant/main'
include { TABIX_TABIX as TABIX_VC_DEEPVARIANT_GVCF  } from '../../../../modules/nf-core/modules/tabix/tabix/main'
include { TABIX_TABIX as TABIX_VC_DEEPVARIANT_VCF   } from '../../../../modules/nf-core/modules/tabix/tabix/main'

//TODO: benchmark if it is better to provide multiple bed files & run on multiple machines + mergeing afterwards || one containing all intervals and run on one larger machine
// Deepvariant: https://github.com/google/deepvariant/issues/510
workflow RUN_DEEPVARIANT {
    take:
    cram                     // channel: [mandatory] [meta, cram, crai, interval]
    dict                     // channel: [optional]
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

    MERGE_DEEPVARIANT_VCF(
        deepvariant_vcf_out.intervals
            .map{ meta, vcf ->

                new_meta = [patient:meta.patient, sample:meta.sample, status:meta.status, gender:meta.gender, id:meta.sample, num_intervals:meta.num_intervals]

                [groupKey(new_meta, meta.num_intervals), vcf]
            }.groupTuple(),
        dict)

    MERGE_DEEPVARIANT_GVCF(
        deepvariant_gvcf_out.intervals
            .map{ meta, vcf ->

                new_meta = [patient:meta.patient, sample:meta.sample, status:meta.status, gender:meta.gender, id:meta.sample, num_intervals:meta.num_intervals]

                [groupKey(new_meta, meta.num_intervals), vcf]
            }.groupTuple(),
        dict)

    // Mix output channels for "no intervals" and "with intervals" results
    deepvariant_gvcf = Channel.empty().mix(
                        MERGE_DEEPVARIANT_GVCF.out.vcf,
                        deepvariant_gvcf_out.no_intervals)
                    .map{ meta, vcf ->
                        [[patient:meta.patient, sample:meta.sample, status:meta.status, gender:meta.gender, id:meta.sample, num_intervals:meta.num_intervals, variantcaller:"Deepvariant", type: "gvcf"], vcf]
                    }
    deepvariant_vcf = Channel.empty().mix(
                        MERGE_DEEPVARIANT_VCF.out.vcf,
                        deepvariant_vcf_out.no_intervals)
                    .map{ meta, vcf ->
                        [[patient:meta.patient, sample:meta.sample, status:meta.status, gender:meta.gender, id:meta.sample, num_intervals:meta.num_intervals, variantcaller:"Deepvariant", type: "vcf"], vcf]
                    }

    ch_versions = ch_versions.mix(MERGE_DEEPVARIANT_GVCF.out.versions)
    ch_versions = ch_versions.mix(MERGE_DEEPVARIANT_VCF.out.versions)
    ch_versions = ch_versions.mix(DEEPVARIANT.out.versions)
    ch_versions = ch_versions.mix(TABIX_VC_DEEPVARIANT_GVCF.out.versions)
    ch_versions = ch_versions.mix(TABIX_VC_DEEPVARIANT_VCF.out.versions)

    emit:
    deepvariant_vcf
    deepvariant_gvcf
    versions = ch_versions
}
