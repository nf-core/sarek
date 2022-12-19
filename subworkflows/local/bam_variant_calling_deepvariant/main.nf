include { GATK4_MERGEVCFS as MERGE_DEEPVARIANT_GVCF } from '../../../modules/nf-core/gatk4/mergevcfs/main'
include { GATK4_MERGEVCFS as MERGE_DEEPVARIANT_VCF  } from '../../../modules/nf-core/gatk4/mergevcfs/main'
include { DEEPVARIANT                               } from '../../../modules/nf-core/deepvariant/main'
include { TABIX_TABIX as TABIX_VC_DEEPVARIANT_GVCF  } from '../../../modules/nf-core/tabix/tabix/main'
include { TABIX_TABIX as TABIX_VC_DEEPVARIANT_VCF   } from '../../../modules/nf-core/tabix/tabix/main'

// Deepvariant: https://github.com/google/deepvariant/issues/510
workflow BAM_VARIANT_CALLING_DEEPVARIANT {
    take:
    cram                     // channel: [mandatory] [meta, cram, crai, interval]
    dict                     // channel: [optional]
    fasta                    // channel: [mandatory]
    fasta_fai                // channel: [mandatory]

    main:
    versions = Channel.empty()

    DEEPVARIANT(cram, fasta, fasta_fai)

    vcf_out = DEEPVARIANT.out.vcf.branch{
        intervals:    it[0].num_intervals > 1
        no_intervals: it[0].num_intervals <= 1
    }

    gvcf_out = DEEPVARIANT.out.gvcf.branch{
        intervals:    it[0].num_intervals > 1
        no_intervals: it[0].num_intervals <= 1
    }

    // Only when no_intervals
    TABIX_VC_DEEPVARIANT_VCF(vcf_out.no_intervals)
    TABIX_VC_DEEPVARIANT_GVCF(gvcf_out.no_intervals)

    // Only when using intervals

    MERGE_DEEPVARIANT_VCF(
        vcf_out.intervals.map{ meta, vcf ->
            [ groupKey(meta.subMap('num_intervals', 'patient', 'sample', 'sex', 'status')
                + [ id: meta.sample ], meta.num_intervals), vcf ]
        }.groupTuple(),
        dict.map{ it -> [[id:it[0].baseName], it]})

    MERGE_DEEPVARIANT_GVCF(
        gvcf_out.intervals.map{ meta, vcf ->
            [ groupKey(meta.subMap('num_intervals', 'patient', 'sample', 'sex', 'status')
                + [ id: meta.sample ], meta.num_intervals), vcf ]
        }.groupTuple(),
        dict.map{ it -> [[id:it[0].baseName], it]})

    // Mix output channels for intervals and no_intervals results
    gvcf = Channel.empty().mix(MERGE_DEEPVARIANT_GVCF.out.vcf, gvcf_out.no_intervals).map{ meta, vcf ->
        [ meta.subMap('num_intervals', 'patient', 'sample', 'sex', 'status')
            + [ id: meta.sample, variantcaller:"deepvariant" ],
        vcf ]
    }

    // Mix output channels for intervals and no_intervals results
    vcf = Channel.empty().mix(MERGE_DEEPVARIANT_VCF.out.vcf, vcf_out.no_intervals).map{ meta, vcf ->
        [ meta.subMap('num_intervals', 'patient', 'sample', 'sex', 'status')
            + [ id: meta.sample, variantcaller:"deepvariant" ],
        vcf ]
    }

    versions = versions.mix(MERGE_DEEPVARIANT_GVCF.out.versions)
    versions = versions.mix(MERGE_DEEPVARIANT_VCF.out.versions)
    versions = versions.mix(DEEPVARIANT.out.versions)
    versions = versions.mix(TABIX_VC_DEEPVARIANT_GVCF.out.versions)
    versions = versions.mix(TABIX_VC_DEEPVARIANT_VCF.out.versions)

    emit:
    gvcf
    vcf

    versions
}
