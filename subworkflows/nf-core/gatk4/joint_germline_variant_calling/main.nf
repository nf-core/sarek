//
// merge samples with genomicsdbimport, perform joint genotyping with genotypeGVCFS
include { GATK4_GENOMICSDBIMPORT }                 from '../../../../modules/nf-core/modules/gatk4/genomicsdbimport/main'
include { GATK4_GENOTYPEGVCFS }                    from '../../../../modules/nf-core/modules/gatk4/genotypegvcfs/main'
include { GATK4_MERGEVCFS as MERGE_GENOTYPEGVCFS } from '../../../../modules/nf-core/modules/gatk4/mergevcfs/main'
include { GATK4_VARIANTRECALIBRATOR as SNP_VQSR }  from '../../../../modules/nf-core/modules/gatk4/variantrecalibrator/main'
include { GATK4_VARIANTRECALIBRATOR as INDEL_VQSR} from '../../../../modules/nf-core/modules/gatk4/variantrecalibrator/main'

workflow GATK_JOINT_GERMLINE_VARIANT_CALLING {
    take:
    input            // channel: [ val(meta), [ input ], [ input_index ], interval, [], []]
    fasta            // channel: /path/to/reference/fasta
    fai              // channel: /path/to/reference/fasta/index
    dict             // channel: /path/to/reference/fasta/dictionary
    sites            // channel: /path/to/known/sites/file
    sites_index      // channel: /path/to/known/sites/index

    main:
    ch_versions    = Channel.empty()
    resources_snp   = (params.RESOURCE_SNP     ? params.genome.RESOURCE_SNP   : Channel.empty() )
    resources_indel = (params.RESOURCE_INDEL   ? params.genome.RESOURCE_INDEL : Channel.empty() )

    gendb_input = input.map{
        meta, gvcf, tbi ->
            interval_file = meta.num_intervals > 1 ? []      : params.intervals
            interval_val  = meta.num_intervals > 1 ? meta.id : []
            [meta, gvcf, tbi, interval_file, interval_val, []]
        }

    //
    //Convert all sample vcfs into a genomicsdb workspace using genomicsdbimport.
    //
    GATK4_GENOMICSDBIMPORT ( gendb_input, false, false, false )
    genotype_input = GATK4_GENOMICSDBIMPORT.out.genomicsdb.join(gendb_input).map{
        meta, genomicsdb, gvcf, tbi, interval_file, interval_val, wpath ->
            new_meta = meta
            new_meta.id = meta.id.replaceAll(":","_")
            [new_meta, genomicsdb, [], [], []]
        }
    ch_versions = ch_versions.mix(GATK4_GENOMICSDBIMPORT.out.versions)

    //
    //Joint genotyping performed using GenotypeGVCFs
    //
    vcfs = GATK4_GENOTYPEGVCFS ( genotype_input, fasta, fai, dict, sites, sites_index).vcf
    merge_vcfs_input = vcfs.map { meta, vcf ->
        [[id:"joint variant calling", num_intervals: meta.num_intervals], vcf]}.groupTuple()

    //
    //Merge vcfs called by interval into a single VCF
    //
    merged_vcf = MERGE_GENOTYPEGVCFS(merge_vcfs_input,dict)

    vqsr_input = merged_vcf.vcf.join(merged_vcf.tbi)
    SNP_VQSR(vqsr_input,
             resources_snp,
             fasta,
             fai,
             dict)
    INDEL_VQSR(vqsr_input,
             resources_indel,
             fasta,
             fai,
             dict)

    ch_versions = ch_versions.mix(GATK4_GENOTYPEGVCFS.out.versions)

    emit:
    versions       = ch_versions                           // channel: [ versions.yml ]
    genomicsdb     = GATK4_GENOMICSDBIMPORT.out.genomicsdb // channel: [ val(meta), [ genomicsdb ] ]
    genotype_vcf   = MERGE_GENOTYPEGVCFS.out.vcf           // channel: [ val(meta), [ vcf ] ]
    genotype_index = MERGE_GENOTYPEGVCFS.out.tbi           // channel: [ val(meta), [ tbi ] ]
}
