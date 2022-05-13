//
// merge samples with genomicsdbimport, perform joint genotyping with genotypeGVCFS


include { GATK4_GENOMICSDBIMPORT } from '../../../modules/gatk4/genomicsdbimport/main'
include { GATK4_GENOTYPEGVCFS } from '../../../modules/gatk4/genotypegvcfs/main'

workflow GATK_JOINT_GERMLINE_VARIANT_CALLING {
    take:
    input            // channel: [ val(meta), [ input ], [ input_index ], [] ]
    fasta            // channel: /path/to/reference/fasta
    fai              // channel: /path/to/reference/fasta/index
    dict             // channel: /path/to/reference/fasta/dictionary
    sites            // channel: /path/to/known/sites/file
    sites_index      // channel: /path/to/known/sites/index

    main:
    ch_versions = Channel.empty()
    ch_joint_intervals = input.map{ meta, vcf, tbi, intervals, interval, dummy -> [meta, intervals, interval]  }

    //
    //Convert all sample vcfs into a genomicsdb workspace using genomicsdbimport.
    //
    gendb_input = input
    GATK4_GENOMICSDBIMPORT ( gendb_input, false, false, false )

    ch_versions = ch_versions.mix(GATK4_GENOMICSDBIMPORT.out.versions)

    //
    //Joint genotyping performed using GenotypeGVCFs
    //
    //[] is a placeholder for the input where the vcf tbi file would be passed in for non-genomicsdb workspace runs, which is not necessary for this workflow as it uses genomicsdb workspaces.
    genotype_input = GATK4_GENOMICSDBIMPORT.out.genomicsdb.combine(ch_joint_intervals, by: 0)
    GATK4_GENOTYPEGVCFS ( genotype_input, fasta, fai, dict, sites, sites_index)

    ch_versions = ch_versions.mix(GATK4_GENOTYPEGVCFS.out.versions)

    emit:
    versions       = ch_versions                           // channel: [ versions.yml ]
    genomicsdb     = GATK4_GENOMICSDBIMPORT.out.genomicsdb // channel: [ val(meta), [ genomicsdb ] ]
    genotype_vcf   = GATK4_GENOTYPEGVCFS.out.vcf           // channel: [ val(meta), [ vcf ] ]
    genotype_index = GATK4_GENOTYPEGVCFS.out.tbi           // channel: [ val(meta), [ tbi ] ]
}
