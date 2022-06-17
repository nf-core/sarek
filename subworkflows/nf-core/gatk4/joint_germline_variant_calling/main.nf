//
// merge samples with genomicsdbimport, perform joint genotyping with genotypeGVCFS
include { GATK4_GENOMICSDBIMPORT } from '../../../../modules/nf-core/modules/gatk4/genomicsdbimport/main'
include { GATK4_GENOTYPEGVCFS } from '../../../../modules/nf-core/modules/gatk4/genotypegvcfs/main'
include { GATK4_MERGEVCFS as MERGE_GENOTYPEGVCFS } from '../../../../modules/nf-core/modules/gatk4/mergevcfs/main'

workflow GATK_JOINT_GERMLINE_VARIANT_CALLING {
    take:
    input            // channel: [ val(meta), [ input ], [ input_index ], interval, [], []]
    fasta            // channel: /path/to/reference/fasta
    fai              // channel: /path/to/reference/fasta/index
    dict             // channel: /path/to/reference/fasta/dictionary
    sites            // channel: /path/to/known/sites/file
    sites_index      // channel: /path/to/known/sites/index

    main:
    ch_versions = Channel.empty()

    gendb_input = input.map{
        meta, gvcf, tbi ->
            interval_file = meta.num_intervals > 1 ? []                                                    : params.intervals
            interval_val  = meta.num_intervals > 1 ? (gvcf.simpleName - "${meta.id}_").replaceFirst("_",":") : []
            interval_name = meta.num_intervals > 1 ? gvcf.simpleName: meta.id
            new_meta = [patient:meta.patient,
                        sample:meta.sample,
                        id:meta.id,
                        interval_name:interval_name, 
                        num_intervals:meta.num_intervals]
            [new_meta, gvcf, tbi, interval_file, interval_val, []]
        }
    gendb_input.dump(tag:"gdb")

    //
    //Convert all sample vcfs into a genomicsdb workspace using genomicsdbimport.
    //
    GATK4_GENOMICSDBIMPORT ( gendb_input, false, false, false )
    genotype_input = GATK4_GENOMICSDBIMPORT.out.genomicsdb.join(gendb_input).map{
        new_meta, genomicsdb, gvcf, tbi, interval_file, interval_val, wpath ->
            [new_meta, genomicsdb, [], interval_file, interval_val]
        }
    genotype_input.view()
    ch_versions = ch_versions.mix(GATK4_GENOMICSDBIMPORT.out.versions)

    //
    //Joint genotyping performed using GenotypeGVCFs
    //
    vcfs = GATK4_GENOTYPEGVCFS ( genotype_input, fasta, fai, dict, sites, sites_index).vcf
    merge_vcfs_input = vcfs.map{ meta, vcf  ->
        // remove the coordinates from meta to group by patient, sample, ... for merging
        new_meta = [patient:meta.patient, sample:meta.sample, status:meta.status,
                    gender:meta.gender, id:meta.sample, num_intervals:meta.num_intervals]
        [groupKey(new_meta, new_meta.num_intervals), vcf]
    }.groupTuple()

   //
   //Merge vcfs called by interval into a single VCF
   //
   MERGE_GENOTYPEGVCFS(merge_vcfs_input,dict)


    ch_versions = ch_versions.mix(GATK4_GENOTYPEGVCFS.out.versions)

    emit:
    versions       = ch_versions                           // channel: [ versions.yml ]
    genomicsdb     = GATK4_GENOMICSDBIMPORT.out.genomicsdb // channel: [ val(meta), [ genomicsdb ] ]
    genotype_vcf   = MERGE_GENOTYPEGVCFS.out.vcf           // channel: [ val(meta), [ vcf ] ]
    genotype_index = MERGE_GENOTYPEGVCFS.out.tbi           // channel: [ val(meta), [ tbi ] ]
}
