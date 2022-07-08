//
// merge samples with genomicsdbimport, perform joint genotyping with genotypeGVCFS
include { BCFTOOLS_SORT }                          from '../../../../modules/nf-core/modules/bcftools/sort/main'
include { GATK4_GENOMICSDBIMPORT }                 from '../../../../modules/nf-core/modules/gatk4/genomicsdbimport/main'
include { GATK4_GENOTYPEGVCFS }                    from '../../../../modules/nf-core/modules/gatk4/genotypegvcfs/main'
include { GATK4_MERGEVCFS as MERGE_GENOTYPEGVCFS } from '../../../../modules/nf-core/modules/gatk4/mergevcfs/main'
include { GATK4_MERGEVCFS as MERGE_VQSR }          from '../../../../modules/nf-core/modules/gatk4/mergevcfs/main'
include { GATK4_VARIANTRECALIBRATOR as VARIANTRECALIBRATOR_SNP }  from '../../../../modules/nf-core/modules/gatk4/variantrecalibrator/main'
include { GATK4_APPLYVQSR as GATK4_APPLYVQSR_SNP } from '../../../../modules/nf-core/modules/gatk4/applyvqsr/main'
include { GATK4_APPLYVQSR as GATK4_APPLYVQSR_INDEL } from '../../../../modules/nf-core/modules/gatk4/applyvqsr/main'
include { GATK4_VARIANTRECALIBRATOR as VARIANTRECALIBRATOR_INDEL} from '../../../../modules/nf-core/modules/gatk4/variantrecalibrator/main'

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
    resources_snps = Channel.empty()
    resources_indels = Channel.empty()


    known_snps      = params.known_snps      ? Channel.fromPath(params.known_snps).collect() : Channel.empty()
    known_snps_tbi  = params.known_snps_tbi  ? Channel.fromPath(params.known_snps_tbi).collect() : Channel.empty()


    known_indels       = params.known_indels ? Channel.fromPath(params.known_indels).collect() : Channel.empty()
    known_indels_tbi   = params.known_indels_tbi ? Channel.fromPath(params.known_indels_tbi).collect() : Channel.empty()

    dbsnp       = params.dbsnp      ? Channel.fromPath(params.dbsnp).collect() : Channel.empty()
    dbsnp_tbi   = params.dbsnp_tbi  ? Channel.fromPath(params.dbsnp_tbi).collect() : Channel.empty()

    resource_snps_vcf = known_snps.combine(dbsnp)
    resource_snps_tbi = known_snps_tbi.combine(dbsnp_tbi)

    resource_indels_vcf = known_indels.combine(dbsnp)
    resource_indels_tbi = known_indels_tbi.combine(dbsnp_tbi)
    
    gendb_input = input.map{
        meta, gvcf, tbi, intervals->
            new_meta = [id: meta.num_intervals > 1 ? "joint variant calling" : "no_intervals",
                        intervals_name: meta.intervals_name,
                        num_intervals: meta.num_intervals
                       ]
            [ new_meta, gvcf, tbi, intervals ]
        }.groupTuple(by:[0,3]).map{ new_meta, gvcf, tbi, intervals -> 
            interval_file = new_meta.num_intervals > 1 ? intervals : params.intervals
            [ new_meta, gvcf, tbi, interval_file, [], [] ] }

    //
    //Convert all sample vcfs into a genomicsdb workspace using genomicsdbimport.
    //
    GATK4_GENOMICSDBIMPORT ( gendb_input, false, false, false )

    genotype_input = GATK4_GENOMICSDBIMPORT.out.genomicsdb.map{
        meta, genomicsdb ->
            [meta, genomicsdb, [], [], []]
        }
    ch_versions = ch_versions.mix(GATK4_GENOMICSDBIMPORT.out.versions)

    //
    //Joint genotyping performed using GenotypeGVCFs
    //Sort vcfs called by interval within each VCF
    //

    vcfs = GATK4_GENOTYPEGVCFS ( genotype_input, fasta, fai, dict, sites, sites_index).vcf

    vcfs_sorted_input = BCFTOOLS_SORT(vcfs).vcf
    merge_vcfs_sorted_input = vcfs_sorted_input.map { meta, vcf ->
        [[id:"joint_variant_calling", num_intervals: meta.num_intervals], vcf]}.groupTuple()

    //
    //Merge vcfs called by interval into a single VCF
    //
    merged_vcf = MERGE_GENOTYPEGVCFS(merge_vcfs_sorted_input,dict)

    if (params.variant_recalibration){
        vqsr_input = merged_vcf.vcf.join(merged_vcf.tbi)

        //
        //Recalibrate SNP and INDEL separately.
        //
        VARIANTRECALIBRATOR_SNP(vqsr_input,
                resource_snps_vcf,
                resource_snps_tbi,
                fasta,
                fai,
                dict)

        VARIANTRECALIBRATOR_INDEL(vqsr_input,
                resource_indels_vcf,
                resource_indels_tbi,
                fasta,
                fai,
                dict)

        //
        //Prepare SNP and INDEL separately for ApplyVQSR
        //
        vqsr_input.dump(tag:"in")
        VARIANTRECALIBRATOR_INDEL.out.recal.dump(tag:"vqsr")
       
        vqsr_input_snp   = vqsr_input.join( VARIANTRECALIBRATOR_SNP.out.recal).join(
                                            VARIANTRECALIBRATOR_SNP.out.idx).join(
                                            VARIANTRECALIBRATOR_SNP.out.tranches)

        vqsr_input_indel   = vqsr_input.join( VARIANTRECALIBRATOR_INDEL.out.recal).join(
                                            VARIANTRECALIBRATOR_INDEL.out.idx).join(
                                            VARIANTRECALIBRATOR_INDEL.out.tranches)

        vqsr_snp_vcf = GATK4_APPLYVQSR_SNP(vqsr_input_snp,
                            fasta,
                            fai,
                            dict ).vcf
                         
        vqsr_indel_vcf = GATK4_APPLYVQSR_INDEL(vqsr_input_indel,
                            fasta,
                            fai,
                            dict ).vcf

        vqsr_snp_vcf.mix(vqsr_indel_vcf).groupTuple().dump(tag:"merge")
        //
        //Merge VQSR outputs into final VCF
        //
        output_ch = MERGE_VQSR(
            vqsr_snp_vcf.mix(vqsr_indel_vcf).groupTuple(),
            dict
        )

    } else {
        output_ch = merged_vcf
    }

    ch_versions = ch_versions.mix(GATK4_GENOTYPEGVCFS.out.versions,
                                  VARIANTRECALIBRATION_SNP.versions,
                                  GATK4_APPLYVQSR_SNP.out.versions
                                 )

    emit:
    versions       = ch_versions                           // channel: [ versions.yml ]
    genotype_vcf   = output_ch.vcf           // channel: [ val(meta), [ vcf ] ]
    genotype_index = output_ch.tbi           // channel: [ val(meta), [ tbi ] ]
}
