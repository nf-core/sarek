include { GATK4_MERGEVCFS as MERGE_HAPLOTYPECALLER } from '../../../../modules/nf-core/modules/gatk4/mergevcfs/main'
include { GATK4_GENOTYPEGVCFS as GENOTYPEGVCFS     } from '../../../../modules/nf-core/modules/gatk4/genotypegvcfs/main'
include { GATK4_HAPLOTYPECALLER as HAPLOTYPECALLER } from '../../../../modules/nf-core/modules/gatk4/haplotypecaller/main'
include { GATK_JOINT_GERMLINE_VARIANT_CALLING      } from '../../../../subworkflows/nf-core/gatk4/joint_germline_variant_calling/main'

workflow RUN_HAPLOTYPECALLER {
    take:
    cram                            // channel: [mandatory] [meta, cram, crai, interval.bed.gz, interval.bed.gz.tbi]
    fasta                           // channel: [mandatory]
    fasta_fai                       // channel: [mandatory]
    dict                            // channel: [mandatory]
    dbsnp                           // channel: [mandatory]
    dbsnp_tbi                       // channel: [mandatory]

    main:

    ch_versions = Channel.empty()

    HAPLOTYPECALLER(
        cram,
        fasta,
        fasta_fai,
        dict,
        dbsnp,
        dbsnp_tbi)



    // Merge after processing
    //single haplotypecaller (default mode)-> CNNScoreVariants
    //joint workflow haplotypecaller (ERC mode) -> GenomicsDBimport -> GenotypeGVCFs -> VQSR

    if (params.joint_germline) {
        genotype_gvcf_to_call = HAPLOTYPECALLER.out.vcf.join(HAPLOTYPECALLER.out.tbi)
        genotype_gvcf_to_call.dump(tag:"htc")

        genotype_vcf = GATK_JOINT_GERMLINE_VARIANT_CALLING(
             genotype_gvcf_to_call,
             fasta,
             fasta_fai,
             dict,
             dbsnp,
             dbsnp_tbi).genotype_vcf

        ch_versions = ch_versions.mix(GATK_JOINT_GERMLINE_VARIANT_CALLING.out.versions)
    } else {
        // Figure out if using intervals or no_intervals
        HAPLOTYPECALLER.out.vcf.branch{
                intervals:    it[0].num_intervals > 1
                no_intervals: it[0].num_intervals <= 1
            }.set{haplotypecaller_vcf_branch}

        HAPLOTYPECALLER.out.tbi.branch{
                intervals:    it[0].num_intervals > 1
                no_intervals: it[0].num_intervals <= 1
            }.set{haplotypecaller_tbi_branch}

        MERGE_HAPLOTYPECALLER(
            haplotypecaller_vcf_branch.intervals
                .map{ meta, vcf ->
                    new_meta = [patient:meta.patient, sample:meta.sample, status:meta.status, gender:meta.gender, id:meta.sample, num_intervals:meta.num_intervals]
                [groupKey(new_meta, new_meta.num_intervals), vcf]
            }.groupTuple(),
            dict)

        haplotypecaller_vcf = Channel.empty().mix(
            MERGE_HAPLOTYPECALLER.out.vcf,
            haplotypecaller_vcf_branch.no_intervals)

        haplotypecaller_tbi = Channel.empty().mix(
            MERGE_HAPLOTYPECALLER.out.tbi,
            haplotypecaller_vcf_branch.no_intervals)

        // CNNScoreVariants
    }

    genotype_vcf = Channel.empty()
    haplotypecaller_vcf = Channel.empty()
    //ch_versions = ch_versions.mix(MERGE_HAPLOTYPECALLER.out.versions)
    //ch_versions = ch_versions.mix(GENOTYPEGVCFS.out.versions)
    //ch_versions = ch_versions.mix(GATK_JOINT_GERMLINE_VARIANT_CALLING.out.versions)
    ch_versions = ch_versions.mix(HAPLOTYPECALLER.out.versions)
    // ch_versions = ch_versions.mix(TABIX_VC_HAPLOTYPECALLER.out.versions)

    emit:
    versions = ch_versions
    genotype_vcf
    haplotypecaller_vcf
}
