include { BCFTOOLS_SORT                                } from '../../../../modules/nf-core/modules/bcftools/sort/main'
include { GATK4_MERGEVCFS as MERGE_FREEBAYES           } from '../../../../modules/nf-core/modules/gatk4/mergevcfs/main'
include { FREEBAYES                                    } from '../../../../modules/nf-core/modules/freebayes/main'
include { TABIX_TABIX as TABIX_VC_FREEBAYES            } from '../../../../modules/nf-core/modules/tabix/tabix/main'

workflow RUN_FREEBAYES {
    take:
    cram                     // channel: [mandatory] [meta, cram, crai, [], [], interval]
    dict
    fasta                    // channel: [mandatory]
    fasta_fai                // channel: [mandatory]

    main:

    ch_versions = Channel.empty()

    FREEBAYES(
        cram,
        fasta,
        fasta_fai,
        [], [], [])

    BCFTOOLS_SORT(FREEBAYES.out.vcf)
    BCFTOOLS_SORT.out.vcf.branch{
            intervals:    it[0].num_intervals > 1
            no_intervals: it[0].num_intervals <= 1
        }.set{bcftools_vcf_out}

    // Only when no intervals
    TABIX_VC_FREEBAYES(bcftools_vcf_out.no_intervals)

    // Only when using intervals
    MERGE_FREEBAYES(
        bcftools_vcf_out.intervals
            .map{ meta, vcf ->

                new_id = meta.tumor_id ? meta.tumor_id + "_vs_" + meta.normal_id : meta.sample

                new_meta = meta.tumor_id ? [patient:meta.patient, normal_id:meta.normal_id, tumor_id:meta.tumor_id, gender:meta.gender, id:new_id, num_intervals:meta.num_intervals]
                                        : [patient:meta.patient, sample:meta.sample, status:meta.status, gender:meta.gender, id:new_id, num_intervals:meta.num_intervals]
                [groupKey(new_meta, meta.num_intervals), vcf]
            }.groupTuple(),
        dict
    )

    // Mix output channels for "no intervals" and "with intervals" results
    freebayes_vcf = Channel.empty().mix(
                        MERGE_FREEBAYES.out.vcf,
                        bcftools_vcf_out.no_intervals)
                    .map{ meta, vcf ->

                        new_id = meta.tumor_id ? meta.tumor_id + "_vs_" + meta.normal_id : meta.sample

                        new_meta = meta.tumor_id ? [patient:meta.patient, normal_id:meta.normal_id, tumor_id:meta.tumor_id, gender:meta.gender, id:new_id, num_intervals:meta.num_intervals, variantcaller:"Freebayes"]
                                        : [patient:meta.patient, sample:meta.sample, status:meta.status, gender:meta.gender, id:new_id, num_intervals:meta.num_intervals, variantcaller:"Freebayes"]
                        [new_meta, vcf]
                    }

    ch_versions = ch_versions.mix(BCFTOOLS_SORT.out.versions)
    ch_versions = ch_versions.mix(MERGE_FREEBAYES.out.versions)
    ch_versions = ch_versions.mix(FREEBAYES.out.versions)
    ch_versions = ch_versions.mix(TABIX_VC_FREEBAYES.out.versions)

    emit:
    freebayes_vcf
    versions = ch_versions
}
