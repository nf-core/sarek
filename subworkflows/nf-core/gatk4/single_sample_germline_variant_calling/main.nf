include { GATK4_CNNSCOREVARIANTS      as CNNSCOREVARIANTS      } from '../../../../modules/nf-core/modules/gatk4/cnnscorevariants/main'
include { GATK4_FILTERVARIANTTRANCHES as FILTERVARIANTTRANCHES } from  '../../../../modules/nf-core/modules/gatk4/filtervarianttranches/main'
include { GATK4_MERGEVCFS             as MERGE_HAPLOTYPECALLER_FILTERED        } from '../../../../modules/nf-core/modules/gatk4/mergevcfs/main'

workflow GATK_SINGLE_SAMPLE_GERMLINE_VARIANT_CALLING{

    take:
    vcf // meta, vcf, tbi, intervals
    fasta
    fasta_fai
    dict
    known_sites
    known_sites_tbi

    main:

    ch_versions = Channel.empty()

    CNNSCOREVARIANTS(vcf.map{ meta, vcf, tbi, intervals -> [meta,vcf,tbi,[],intervals]},
                    fasta,
                    fasta_fai,
                    dict,
                    [],
                    [])

    cnn_out = CNNSCOREVARIANTS.out.vcf.join(CNNSCOREVARIANTS.out.tbi).join(vcf)
        .map{   meta, cnn_vcf,cnn_tbi, haplotyc_vcf, haplotyc_tbi, intervals
            -> [meta, cnn_vcf, cnn_tbi, intervals]
        }

    cnn_out.view()

    // FILTERVARIANTTRANCHES(cnn_out,
    //                         known_sites,
    //                         known_sites_tbi,
    //                         fasta,
    //                         fasta_fai,
    //                         dict)

    //     // Figure out if using intervals or no_intervals
    // HAPLOTYPECALLER.out.vcf.branch{
    //         intervals:    it[0].num_intervals > 1
    //         no_intervals: it[0].num_intervals <= 1
    //     }.set{haplotypecaller_vcf_branch}

    // HAPLOTYPECALLER.out.tbi.branch{
    //         intervals:    it[0].num_intervals > 1
    //         no_intervals: it[0].num_intervals <= 1
    //     }.set{haplotypecaller_tbi_branch}

    // // Only when using intervals
    // MERGE_HAPLOTYPECALLER(
    //     haplotypecaller_vcf_branch.intervals
    //         .map{ meta, vcf ->

    //             new_meta = [patient:meta.patient, sample:meta.sample, status:meta.status, gender:meta.gender, id:meta.sample, num_intervals:meta.num_intervals]

    //             [groupKey(new_meta, new_meta.num_intervals), vcf]
    //         }.groupTuple(),
    //     dict)
    ch_versions = ch_versions.mix(CNNSCOREVARIANTS.out.versions)
    //ch_versions = ch_versions.mix(FILTERVARIANTTRANCHES.out.versions)

    emit:
    versions = ch_versions
    vcf = Channel.empty() //FILTERVARIANTTRANCHES.out.vcf
}
