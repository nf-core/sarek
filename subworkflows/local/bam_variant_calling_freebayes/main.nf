include { BCFTOOLS_SORT                                } from '../../../modules/nf-core/bcftools/sort/main'
include { GATK4_MERGEVCFS as MERGE_FREEBAYES           } from '../../../modules/nf-core/gatk4/mergevcfs/main'
include { FREEBAYES                                    } from '../../../modules/nf-core/freebayes/main'
include { TABIX_TABIX as TABIX_VC_FREEBAYES            } from '../../../modules/nf-core/tabix/tabix/main'

workflow BAM_VARIANT_CALLING_FREEBAYES {
    take:
    cram                     // channel: [mandatory] [meta, cram, crai, [], [], interval]
    dict
    fasta                    // channel: [mandatory]
    fasta_fai                // channel: [mandatory]

    main:
    versions = Channel.empty()

    FREEBAYES(
        cram,
        fasta,
        fasta_fai,
        [], [], [])

    BCFTOOLS_SORT(FREEBAYES.out.vcf)
    bcftools_vcf_out = BCFTOOLS_SORT.out.vcf.branch{
        intervals:    it[0].num_intervals > 1
        no_intervals: it[0].num_intervals <= 1
    }

    // Only when no intervals
    TABIX_VC_FREEBAYES(bcftools_vcf_out.no_intervals)

    // Only when using intervals
    MERGE_FREEBAYES(
        bcftools_vcf_out.intervals
            .map{ meta, vcf ->
                [ groupKey(
                    ( meta.tumor_id ?
                        [ meta.subMap('normal_id', 'num_intervals', 'patient', 'sex', 'tumor_id') + [ id:meta.tumor_id + "_vs_" + meta.normal_id ] ] :
                        [ meta.subMap('num_intervals', 'patient', 'sex', 'status') + [ id:meta.sample ] ] ),
                meta.num_intervals), vcf ]
            }.groupTuple(),
        dict.map{ it -> [ [ id:it[0].baseName ], it] })

    // Mix output channels for "no intervals" and "with intervals" results
    vcf = Channel.empty().mix(
        MERGE_FREEBAYES.out.vcf,
        bcftools_vcf_out.no_intervals).map{ meta, vcf ->
            [ meta.subMap('id', 'normal_id', 'num_intervals', 'patient', 'sex', 'tumor_id') + [ variantcaller:meta."freebayes" ], vcf ]
        }

    versions = versions.mix(BCFTOOLS_SORT.out.versions)
    versions = versions.mix(MERGE_FREEBAYES.out.versions)
    versions = versions.mix(FREEBAYES.out.versions)
    versions = versions.mix(TABIX_VC_FREEBAYES.out.versions)

    emit:
    vcf

    versions
}
