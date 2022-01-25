//
// GERMLINE VARIANT CALLING
//

include { BGZIP as BGZIP_DEEPVARIANT_GVCF             } from '../../modules/local/bgzip'
include { BGZIP as BGZIP_DEEPVARIANT_VCF              } from '../../modules/local/bgzip'
include { BGZIP as BGZIP_FREEBAYES                    } from '../../modules/local/bgzip'
include { BGZIP as BGZIP_HAPLOTYPECALLER              } from '../../modules/local/bgzip'
include { BGZIP as BGZIP_MANTA_SMALL_INDELS           } from '../../modules/local/bgzip'
include { BGZIP as BGZIP_MANTA_SV                     } from '../../modules/local/bgzip'
include { BGZIP as BGZIP_MANTA_DIPLOID                } from '../../modules/local/bgzip'
include { BGZIP as BGZIP_STRELKA                      } from '../../modules/local/bgzip'
include { CONCAT_VCF as CONCAT_VCF_DEEPVARIANT        } from '../../modules/local/concat_vcf/main'
include { CONCAT_VCF as CONCAT_GVCF_DEEPVARIANT       } from '../../modules/local/concat_vcf/main'
include { CONCAT_VCF as CONCAT_VCF_FREEBAYES          } from '../../modules/local/concat_vcf/main'
include { CONCAT_VCF as CONCAT_VCF_HAPLOTYPECALLER    } from '../../modules/local/concat_vcf/main'
include { CONCAT_VCF as CONCAT_VCF_MANTA_SMALL_INDELS } from '../../modules/local/concat_vcf/main'
include { CONCAT_VCF as CONCAT_VCF_MANTA_SV           } from '../../modules/local/concat_vcf/main'
include { CONCAT_VCF as CONCAT_VCF_MANTA_DIPLOID      } from '../../modules/local/concat_vcf/main'
include { CONCAT_VCF as CONCAT_VCF_STRELKA            } from '../../modules/local/concat_vcf/main'
include { DEEPVARIANT                                 } from '../../modules/nf-core/modules/deepvariant/main'
include { FREEBAYES                                   } from '../../modules/nf-core/modules/freebayes/main'
include { GATK_JOINT_GERMLINE_VARIANT_CALLING         } from '../../subworkflows/nf-core/joint_germline_variant_calling/main'
include { GATK4_HAPLOTYPECALLER as HAPLOTYPECALLER    } from '../../modules/nf-core/modules/gatk4/haplotypecaller/main'
include { MANTA_GERMLINE                              } from '../../modules/nf-core/modules/manta/germline/main'
include { STRELKA_GERMLINE                            } from '../../modules/nf-core/modules/strelka/germline/main'
include { TIDDIT_SV                                   } from '../../modules/nf-core/modules/tiddit/sv/main'
include { TABIX_BGZIPTABIX as TABIX_BGZIP_TIDDIT_SV   } from '../../modules/nf-core/modules/tabix/bgziptabix/main'

workflow GERMLINE_VARIANT_CALLING {
    take:
        tools             // Mandatory, list of tools to apply
        cram_recalibrated // channel: [mandatory] cram
        dbsnp             // channel: [mandatory] dbsnp
        dbsnp_tbi         // channel: [mandatory] dbsnp_tbi
        dict              // channel: [mandatory] dict
        fasta             // channel: [mandatory] fasta
        fasta_fai         // channel: [mandatory] fasta_fai
        intervals         // channel: [mandatory] intervals/target regions
        intervals_bed_gz_tbi // channel: [mandatory] intervals/target regions index zipped and indexed
        num_intervals     // val: number of intervals that are used to parallelize exection, either based on capture kit or GATK recommended for WGS
        joint_germline    // val: true/false on whether to run joint_germline calling, only works in combination with haplotypecaller at the moment
        intervals_bed_combine_gz_tbi        // channel: [mandatory] intervals/target regions index zipped and indexed in one file
        //target_bed_gz_tbi // channel: [optional]  target_bed_gz_tbi

    main:

    ch_versions = Channel.empty()

    deepvariant_vcf_gz_tbi      = Channel.empty()
    deepvariant_gvcf_gz_tbi     = Channel.empty()
    freebayes_vcf_gz_tbi        = Channel.empty()
    haplotypecaller_gvcf_gz_tbi = Channel.empty()
    strelka_vcf_gz_tbi          = Channel.empty()

    cram_recalibrated.combine(intervals)
        .map{ meta, cram, crai, intervals ->
            new_meta = meta.clone()
            new_meta.id = meta.sample + "_" + intervals.simpleName
            [new_meta, cram, crai, intervals]
        }.set{cram_recalibrated_intervals}

    cram_recalibrated.combine(intervals_bed_gz_tbi)
        .map{ meta, cram, crai, bed, tbi ->
            new_meta = meta.clone()
            new_meta.id = meta.sample + "_" + bed.simpleName
            [new_meta, cram, crai, bed, tbi]
        }.set{cram_recalibrated_intervals_gz_tbi}

    intervals_bed_combine_gz_tbi.map{ bed_gz, tbi -> [bed_gz]}.set{intervals_bed_combine_gz}

    //TODO: benchmark if it is better to provide multiple bed files & run on multiple machines + mergeing afterwards || one containing all intervals and run on one larger machine

    if (tools.contains('deepvariant')) {
        //TODO: research if multiple targets can be provided: open issue, waiting on answer from maintainers
        // Answer: numshards runs in parallel but on one machine
        DEEPVARIANT(
            cram_recalibrated_intervals,
            fasta,
            fasta_fai)

        BGZIP_DEEPVARIANT_VCF(DEEPVARIANT.out.vcf)
        BGZIP_DEEPVARIANT_GVCF(DEEPVARIANT.out.vcf)

        deepvariant_vcf_to_concat = BGZIP_DEEPVARIANT_VCF.out.vcf.groupTuple(size: num_intervals)
        deepvariant_gvcf_to_concat = BGZIP_DEEPVARIANT_GVCF.out.vcf.groupTuple(size: num_intervals)

        CONCAT_VCF_DEEPVARIANT(deepvariant_vcf_to_concat,fasta_fai, intervals_bed_combine_gz)
        CONCAT_GVCF_DEEPVARIANT(deepvariant_gvcf_to_concat,fasta_fai, intervals_bed_combine_gz)

        deepvariant_vcf_gz_tbi = CONCAT_VCF_DEEPVARIANT.out.vcf
        deepvariant_gvcf_gz_tbi = CONCAT_GVCF_DEEPVARIANT.out.vcf

        ch_versions = ch_versions.mix(DEEPVARIANT.out.versions)
        ch_versions = ch_versions.mix(BGZIP_DEEPVARIANT_VCF.out.versions)
        ch_versions = ch_versions.mix(CONCAT_VCF_DEEPVARIANT.out.versions)
    }

    if (tools.contains('freebayes')){

        cram_recalibrated.combine(intervals).map{ meta, cram, crai, intervals ->
            new_meta = meta.clone()
            new_meta.id = meta.sample + "_" + intervals.simpleName
            [new_meta, cram, crai, [], [], intervals]
        }.set{cram_recalibrated_intervals_freebayes}

        FREEBAYES(
            cram_recalibrated_intervals_freebayes,
            fasta,
            fasta_fai,
            [],
            [],
            []
        )

        BGZIP_FREEBAYES(FREEBAYES.out.vcf)
        freebayes_vcf_to_concat = BGZIP_FREEBAYES.out.vcf.groupTuple(size: num_intervals)

        CONCAT_VCF_FREEBAYES(freebayes_vcf_to_concat,fasta_fai, intervals_bed_combine_gz)
        freebayes_vcf_gz_tbi = CONCAT_VCF_FREEBAYES.out.vcf

        ch_versions = ch_versions.mix(FREEBAYES.out.versions)
        ch_versions = ch_versions.mix(BGZIP_FREEBAYES.out.versions)
        ch_versions = ch_versions.mix(CONCAT_VCF_FREEBAYES.out.versions)
    }

    if (tools.contains('haplotypecaller')) {

        HAPLOTYPECALLER(
                cram_recalibrated_intervals,
                fasta,
                fasta_fai,
                dict,
                dbsnp,
                dbsnp_tbi
        )

        BGZIP_HAPLOTYPECALLER(HAPLOTYPECALLER.out.vcf)

        haplotypecaller_gvcf_to_concat = BGZIP_HAPLOTYPECALLER.out.vcf.groupTuple(size: num_intervals)
        CONCAT_VCF_HAPLOTYPECALLER(haplotypecaller_gvcf_to_concat, fasta_fai, intervals_bed_combine_gz)
        haplotypecaller_gvcf_gz_tbi = CONCAT_VCF_HAPLOTYPECALLER.out.vcf

        ch_versions = ch_versions.mix(HAPLOTYPECALLER.out.versions)
        ch_versions = ch_versions.mix(BGZIP_HAPLOTYPECALLER.out.versions)
        ch_versions = ch_versions.mix(CONCAT_VCF_HAPLOTYPECALLER.out.versions)

        if(joint_germline){
            run_haplotypecaller = false
            run_vqsr            = true //parameter?
            //Waiting on some feedback from gavin
            // GATK_JOINT_GERMLINE_VARIANT_CALLING(
            //     haplotypecaller_vcf_gz_tbi,
            //     run_haplotypecaller,
            //     run_vqsr,
            //     fasta,
            //     fasta_fai,
            //     dict,
            //     dbsnp,
            //     dbsnp_tbi,
            //     "joined",
            //      allelespecific?
            //      resources?
            //      annotation?
            //     "BOTH",
            //     true,
            //     truthsensitivity -> parameter or module?
            // )
            // ch_versions = ch_versions.mix(GATK_JOINT_GERMLINE_VARIANT_CALLING.out.versions)
        }

    }

    if (tools.contains('manta')){
        //TODO: Research if splitting by intervals is ok
        //TODO: merge parallelized vcfs, index them all

        MANTA_GERMLINE(
            cram_recalibrated_intervals_gz_tbi,
            fasta,
            fasta_fai
        )

        manta_candidate_small_indels_vcf_tbi = MANTA_GERMLINE.out.candidate_small_indels_vcf.join(MANTA_GERMLINE.out.candidate_small_indels_vcf_tbi)
        manta_candidate_sv_vcf_tbi           = MANTA_GERMLINE.out.candidate_sv_vcf.join(MANTA_GERMLINE.out.candidate_sv_vcf_tbi)
        manta_diploid_sv_vcf_tbi             = MANTA_GERMLINE.out.diploid_sv_vcf.join(MANTA_GERMLINE.out.diploid_sv_vcf)

        ch_versions = ch_versions.mix(MANTA_GERMLINE.out.versions)
    }

    if (tools.contains('strelka')) {
        //TODO: research if multiple targets can be provided: waiting for reply
        //TODO: Pass over dbsnp/knwon_indels?

        STRELKA_GERMLINE(
            cram_recalibrated_intervals_gz_tbi,
            fasta,
            fasta_fai
            )

        BGZIP_STRELKA(STRELKA_GERMLINE.out.vcf)
        strelka_vcf_to_concat = BGZIP_STRELKA.out.vcf.groupTuple(size: num_intervals)

        CONCAT_VCF_STRELKA(strelka_vcf_to_concat,fasta_fai, intervals_bed_combine_gz)
        strelka_vcf_gz_tbi = CONCAT_VCF_STRELKA.out.vcf

        ch_versions = ch_versions.mix(STRELKA_GERMLINE.out.versions)
        ch_versions = ch_versions.mix(BGZIP_STRELKA.out.versions)
        ch_versions = ch_versions.mix(CONCAT_VCF_STRELKA.out.versions)
    }

    if (tools.contains('tiddit')){
        //TODO: test data not running -> maybe something of with the container...?
        //TODO: Pass over dbsnp/knwon_indels?

        // TIDDIT_SV(
        //     cram_recalibrated,
        //     fasta,
        //     fasta_fai
        // )

        // TABIX_BGZIP_TIDDIT_SV(TIDDIT_SV.out.vcf)
        // tiddit_vcf_gz_tbi = TABIX_BGZIP_TIDDIT_SV.out.gz_tbi
        // tiddit_ploidy = TIDDIT_SV.out.ploidy
        // tiddit_signals = TIDDIT_SV.out.signals
        //tiddit_wig     = TIDDIT_SV.out.wig
        //tiddit_gc_wig  = TIDDIT_SV.out.gc_wig

        //ch_versions = ch_versions.mix(TABIX_BGZIP_TIDDIT_SV.out.versions)
        //ch_versions = ch_versions.mix(TIDDIT_SV.out.versions)
    }

    emit:
    deepvariant_vcf_gz_tbi
    deepvariant_gvcf_gz_tbi
    freebayes_vcf_gz_tbi
    haplotypecaller_gvcf_gz_tbi
    strelka_vcf_gz_tbi

    versions = ch_versions
}
