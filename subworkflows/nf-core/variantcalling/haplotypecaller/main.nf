include { TABIX_BGZIP as BGZIP_VC_HAPLOTYPECALLER  } from '../../../../modules/nf-core/modules/tabix/bgzip/main'
include { CONCAT_VCF as CONCAT_HAPLOTYPECALLER     } from '../../../../modules/local/concat_vcf/main'
include { GATK4_GENOTYPEGVCFS as GENOTYPEGVCFS     } from '../../../../modules/nf-core/modules/gatk4/genotypegvcfs/main'
include { GATK4_HAPLOTYPECALLER as HAPLOTYPECALLER } from '../../../../modules/nf-core/modules/gatk4/haplotypecaller/main'
include { GATK_JOINT_GERMLINE_VARIANT_CALLING      } from '../../../../subworkflows/nf-core/gatk4/joint_germline_variant_calling/main'
include { GATK4_CNNSCOREVARIANTS } from '../modules/nf-core/modules/gatk4/cnnscorevariants/main’
workflow RUN_HAPLOTYPECALLER {
    take:
    cram                            // channel: [mandatory] [meta, cram, crai, interval.bed.gz, interval.bed.gz.tbi]
    fasta                           // channel: [mandatory]
    fasta_fai                       // channel: [mandatory]
    dict                            // channel: [mandatory]
    dbsnp                           // channel: [mandatory]
    dbsnp_tbi                       // channel: [mandatory]
    intervals_bed_gz                // channel: [optional]  Contains a bed.gz file of all intervals combined provided with the cram input(s). Mandatory if interval files are used.
    intervals_bed_combine_gz_tbi    // channel: [optional]  Contains a [bed.gz, bed.gz.tbi ]file of all intervals combined provided with the cram input(s). Mandatory if interval files are used.

    main:

    ch_versions = Channel.empty()

    HAPLOTYPECALLER(
        cram,
        fasta,
        fasta_fai,
        dict,
        dbsnp,
        dbsnp_tbi)

    // Figure out if using intervals or no_intervals
    HAPLOTYPECALLER.out.vcf.branch{
            intervals:    it[0].num_intervals > 1
            no_intervals: it[0].num_intervals <= 1
        }.set{haplotypecaller_vcf_branch}

    HAPLOTYPECALLER.out.tbi.branch{
            intervals:    it[0].num_intervals > 1
            no_intervals: it[0].num_intervals <= 1
        }.set{haplotypecaller_tbi_branch}

    // Only when using intervals
    BGZIP_VC_HAPLOTYPECALLER(haplotypecaller_vcf_branch.intervals)

    CONCAT_HAPLOTYPECALLER(
        BGZIP_VC_HAPLOTYPECALLER.out.output
            .map{ meta, vcf ->
                new_meta = meta.clone()
                new_meta.id = new_meta.sample

                def groupKey = groupKey(meta, meta.num_intervals)
                [new_meta, vcf]
            }.groupTuple(),
        fasta_fai,
        intervals_bed_gz)

    haplotypecaller_vcf = Channel.empty().mix(
        CONCAT_HAPLOTYPECALLER.out.vcf,
        haplotypecaller_vcf_branch.no_intervals)

    haplotypecaller_tbi = Channel.empty().mix(
        CONCAT_HAPLOTYPECALLER.out.tbi,
        haplotypecaller_vcf_branch.no_intervals)

    // genotype_gvcf_to_call = haplotypecaller_gvcf.join(haplotypecaller_gvcf_tbi)
    //     .combine(intervals_bed_combine_gz_tbi)
    //     .map{
    //         meta, gvcf, gvf_tbi, intervals, intervals_tbi ->
    //         new_intervals = intervals.simpleName != "no_intervals" ? intervals : []
    //         new_intervals_tbi = intervals_tbi.simpleName != "no_intervals" ? intervals_tbi : []
    //         [meta, gvcf, gvf_tbi, new_intervals, new_intervals_tbi]
    //     }

    // GENOTYPEGVCFS

    // GENOTYPEGVCFS(
    //     genotype_gvcf_to_call,
    //     fasta,
    //     fasta_fai,
    //     dict,
    //     dbsnp,
    //     dbsnp_tbi)
    //workflow haplotypecaller (default mode)-> CNNScoreVariants
    //workflow haplotypecaller (ERC mode) -> GenomicsDBimport -> GenotypeGVCFs -> VQSR

    //genotype_gvcf = GENOTYPEGVCFS.out.vcf

    if (params.joint_germline) {

        haplotypecaller_vcf_list = haplotypecaller_vcf.toList()
        haplotypecaller_tbi_list = haplotypecaller_vcf.toList()

        joint_germline_vcf_tbi = [ [id: "joint_germline"],
                                    haplotypecaller_vcf_list,
                                    haplotypecaller_tbi_list ]
        // GATK_JOINT_GERMLINE_VARIANT_CALLING(
        //     joint_germline_vcf_tbi,
        //     fasta,
        //     fasta_fai,
        //     intervals,
        //     dict,
        //     dbsnp,
        //     dbsnp_tbi,
        //     allelespecific?
        //     resources?
        //     annotation?
        //     "BOTH",
        //     true,
        //     truthsensitivity -> parameter or module?
        // )
        // ch_versions = ch_versions.mix(GATK_JOINT_GERMLINE_VARIANT_CALLING.out.versions)
    } else {
        // CNNScoreVariants
    }


    ch_versions = ch_versions.mix(BGZIP_VC_HAPLOTYPECALLER.out.versions)
    ch_versions = ch_versions.mix(CONCAT_HAPLOTYPECALLER.out.versions)
    //ch_versions = ch_versions.mix(GENOTYPEGVCFS.out.versions)
    //ch_versions = ch_versions.mix(GATK_JOINT_GERMLINE_VARIANT_CALLING.out.versions)
    ch_versions = ch_versions.mix(HAPLOTYPECALLER.out.versions)
    // ch_versions = ch_versions.mix(TABIX_VC_HAPLOTYPECALLER.out.versions)

    emit:
    versions = ch_versions
    //genotype_gvcf
    haplotypecaller_vcf
}
