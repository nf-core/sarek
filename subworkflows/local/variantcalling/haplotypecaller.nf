include { BGZIP as BGZIP_VC_HAPLOTYPECALLER         } from '../../../modules/local/bgzip'
include { CONCAT_VCF as CONCAT_HAPLOTYPECALLER      } from '../../../modules/local/concat_vcf/main'
include { GATK4_GENOTYPEGVCFS as GENOTYPEGVCFS      } from '../../../modules/nf-core/modules/gatk4/genotypegvcfs/main'
include { GATK4_HAPLOTYPECALLER as HAPLOTYPECALLER  } from '../../../modules/nf-core/modules/gatk4/haplotypecaller/main'
include { GATK_JOINT_GERMLINE_VARIANT_CALLING       } from '../../../subworkflows/nf-core/gatk4/joint_germline_variant_calling/main'
include { TABIX_TABIX as TABIX_VC_HAPLOTYPECALLER   } from '../../../modules/nf-core/modules/tabix/tabix/main'

workflow HAPLOTYPECALLER {
    take:
    cram_recalibrated_intervals
    fasta
    fasta_fai
    dict
    dbsnp
    dbsnp_tbi
    num_intervals
    intervals_bed_combine_gz
    intervals_bed_combine_gz_tbi

    main:

    ch_versions = Channel.empty()

    HAPLOTYPECALLER(
        cram_recalibrated_intervals,
        fasta,
        fasta_fai,
        dict,
        dbsnp,
        dbsnp_tbi)

    // Only when no intervals
    TABIX_VC_HAPLOTYPECALLER(HAPLOTYPECALLER.out.vcf)

    // Only when using intervals
    BGZIP_VC_HAPLOTYPECALLER(HAPLOTYPECALLER.out.vcf)

    CONCAT_HAPLOTYPECALLER(
        BGZIP_VC_HAPLOTYPECALLER.out.vcf
            .map{ meta, vcf ->
                new_meta = meta.clone()
                new_meta.id = new_meta.sample
                [new_meta, vcf]
            }.groupTuple(size: num_intervals),
        fasta_fai,
        intervals_bed_combine_gz)

    HAPLOTYPECALLER.out.vcf.groupTuple(size: num_intervals)
        .branch{
            intervals:    it[1].size() > 1
            no_intervals: it[1].size() == 1
        }.set{haplotypecaller_gvcf_intervals}

    HAPLOTYPECALLER.out.tbi.groupTuple(size: num_intervals)
        .branch{
            intervals:    it[1].size() > 1
            no_intervals: it[1].size() == 1
        }.set{haplotypecaller_gvcf_tbi_intervals}

    haplotypecaller_gvcf = Channel.empty().mix(
        CONCAT_HAPLOTYPECALLER.out.vcf,
        haplotypecaller_gvcf_intervals.no_intervals)

    haplotypecaller_gvcf_tbi = Channel.empty().mix(
        CONCAT_HAPLOTYPECALLER.out.tbi,
        haplotypecaller_gvcf_tbi_intervals.no_intervals)

    genotype_gvcf_to_call = haplotypecaller_gvcf.join(haplotypecaller_gvcf_tbi)
        .combine(intervals_bed_combine_gz_tbi)
        .map{
            meta, gvcf, gvf_tbi, intervals, intervals_tbi ->
            new_intervals = intervals.simpleName != "no_intervals" ? intervals : []
            new_intervals_tbi = intervals_tbi.simpleName != "no_intervals" ? intervals_tbi : []
            [meta, gvcf, gvf_tbi, new_intervals, new_intervals_tbi]
        }

    // GENOTYPEGVCFS

    GENOTYPEGVCFS(
        genotype_gvcf_to_call,
        fasta,
        fasta_fai,
        dict,
        dbsnp,
        dbsnp_tbi)

    genotype_gvcf = GENOTYPEGVCFS.out.vcf

    // if (joint_germline) {
    //     run_haplotypecaller = false
    //     run_vqsr            = true //parameter?
    //     some feedback from gavin
    //     GATK_JOINT_GERMLINE_VARIANT_CALLING(
    //         haplotypecaller_vcf_gz_tbi,
    //         run_haplotypecaller,
    //         run_vqsr,
    //         fasta,
    //         fasta_fai,
    //         dict,
    //         dbsnp,
    //         dbsnp_tbi,
    //         "joined",
    //          allelespecific?
    //          resources?
    //          annotation?
    //         "BOTH",
    //         true,
    //         truthsensitivity -> parameter or module?
    //     )
    //     ch_versions = ch_versions.mix(GATK_JOINT_GERMLINE_VARIANT_CALLING.out.versions)
    // }
    ch_versions = ch_versions.mix(BGZIP_VC_HAPLOTYPECALLER.out.versions)
    ch_versions = ch_versions.mix(CONCAT_HAPLOTYPECALLER.out.versions)
    ch_versions = ch_versions.mix(GENOTYPEGVCFS.out.versions)
    ch_versions = ch_versions.mix(HAPLOTYPECALLER.out.versions)
    ch_versions = ch_versions.mix(TABIX_VC_HAPLOTYPECALLER.out.versions)

    emit:
    versions = ch_versions
    genotype_gvcf
    haplotypecaller_gvcf
}
