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
include { TABIX_TABIX as TABIX_DEEPVARIANT_VCF        } from '../../modules/nf-core/modules/tabix/tabix/main'
include { TABIX_TABIX as TABIX_DEEPVARIANT_GVCF       } from '../../modules/nf-core/modules/tabix/tabix/main'
include { TABIX_TABIX as TABIX_FREEBAYES              } from '../../modules/nf-core/modules/tabix/tabix/main'
include { TABIX_TABIX as TABIX_HAPLOTYPECALLER        } from '../../modules/nf-core/modules/tabix/tabix/main'
include { TABIX_TABIX as TABIX_MANTA                  } from '../../modules/nf-core/modules/tabix/tabix/main'
include { TABIX_TABIX as TABIX_STRELKA                } from '../../modules/nf-core/modules/tabix/tabix/main'


workflow GERMLINE_VARIANT_CALLING {
    take:
        tools                               // Mandatory, list of tools to apply
        cram_recalibrated                   // channel: [mandatory] cram
        dbsnp                               // channel: [mandatory] dbsnp
        dbsnp_tbi                           // channel: [mandatory] dbsnp_tbi
        dict                                // channel: [mandatory] dict
        fasta                               // channel: [mandatory] fasta
        fasta_fai                           // channel: [mandatory] fasta_fai
        intervals                           // channel: [mandatory] intervals/target regions
        intervals_bed_gz_tbi                // channel: [mandatory] intervals/target regions index zipped and indexed
        intervals_bed_combine_gz_tbi        // channel: [mandatory] intervals/target regions index zipped and indexed in one file
        num_intervals                       // val: number of intervals that are used to parallelize exection, either based on capture kit or GATK recommended for WGS
        no_intervals
        joint_germline                      // val: true/false on whether to run joint_germline calling, only works in combination with haplotypecaller at the moment

    main:

    if(!tools) tools = ""

    ch_versions = Channel.empty()

    deepvariant_vcf_gz_tbi                  = Channel.empty()
    deepvariant_gvcf_gz_tbi                 = Channel.empty()
    freebayes_vcf_gz_tbi                    = Channel.empty()
    haplotypecaller_gvcf_gz_tbi             = Channel.empty()
    manta_candidate_small_indels_vcf_tbi    = Channel.empty()
    manta_candidate_sv_vcf_tbi              = Channel.empty()
    manta_diploid_sv_vcf_tbi                = Channel.empty()
    strelka_vcf_gz_tbi                      = Channel.empty()

    cram_recalibrated.combine(intervals)
        .map{ meta, cram, crai, intervals ->
            new_meta = meta.clone()
            new_meta.id = intervals.baseName != "no_intervals" ? meta.sample + "_" + intervals.baseName : meta.sample
            intervals = intervals.baseName != "no_intervals" ? intervals : []
            [new_meta, cram, crai, intervals]
        }.set{cram_recalibrated_intervals}

    cram_recalibrated_intervals.view()

    cram_recalibrated.combine(intervals_bed_gz_tbi)
        .map{ meta, cram, crai, bed, tbi ->
            new_meta = meta.clone()
            new_meta.id = bed.simpleName != "no_intervals" ? meta.sample + "_" + bed.simpleName : meta.sample
            bed = bed.simpleName != "no_intervals" ? bed : []
            tbi = tbi.simpleName != "no_intervals" ? tbi : []

            [new_meta, cram, crai, bed, tbi]
        }.set{cram_recalibrated_intervals_gz_tbi}

    //TODO PASS THIS OVER (ALREADY PROVIDED IN SARKE..NF)
    intervals_bed_combine_gz_tbi.map{ bed_gz, tbi -> [bed_gz]}.set{intervals_bed_combine_gz}

    //TODO: benchmark if it is better to provide multiple bed files & run on multiple machines + mergeing afterwards || one containing all intervals and run on one larger machine
    // Deepvariant: https://github.com/google/deepvariant/issues/510

    if (tools.contains('deepvariant')) {

        DEEPVARIANT(
            cram_recalibrated_intervals,
            fasta,
            fasta_fai)
        ch_versions = ch_versions.mix(DEEPVARIANT.out.versions)

        if(no_intervals){
            TABIX_DEEPVARIANT_VCF(DEEPVARIANT.out.vcf)
            TABIX_DEEPVARIANT_GVCF(DEEPVARIANT.out.gvcf)

            deepvariant_vcf_gz_tbi = DEEPVARIANT.out.vcf.join(TABIX_DEEPVARIANT_VCF.out.tbi)
            deepvariant_gvcf_gz_tbi = DEEPVARIANT.out.gvcf.join(TABIX_DEEPVARIANT_GVCF.out.tbi)

            ch_versions = ch_versions.mix(TABIX_DEEPVARIANT_VCF.out.versions)
            ch_versions = ch_versions.mix(TABIX_DEEPVARIANT_GVCF.out.versions)
        }else{
            BGZIP_DEEPVARIANT_VCF(DEEPVARIANT.out.vcf)
            BGZIP_DEEPVARIANT_GVCF(DEEPVARIANT.out.gvcf)

            deepvariant_vcf_to_concat = BGZIP_DEEPVARIANT_VCF.out.vcf.groupTuple(size: num_intervals)
            deepvariant_gvcf_to_concat = BGZIP_DEEPVARIANT_GVCF.out.vcf.groupTuple(size: num_intervals)

            CONCAT_VCF_DEEPVARIANT(deepvariant_vcf_to_concat,fasta_fai, intervals_bed_combine_gz)
            CONCAT_GVCF_DEEPVARIANT(deepvariant_gvcf_to_concat,fasta_fai, intervals_bed_combine_gz)

            deepvariant_vcf_gz_tbi = CONCAT_VCF_DEEPVARIANT.out.vcf
            deepvariant_gvcf_gz_tbi = CONCAT_GVCF_DEEPVARIANT.out.vcf

            ch_versions = ch_versions.mix(BGZIP_DEEPVARIANT_VCF.out.versions)
            ch_versions = ch_versions.mix(CONCAT_VCF_DEEPVARIANT.out.versions)
        }

    }

    if (tools.contains('freebayes')){

        cram_recalibrated.combine(intervals).map{ meta, cram, crai, intervals ->
            new_meta = meta.clone()
            new_meta.id = meta.sample + "_" + intervals.simpleName
            new_meta.id = intervals.baseName != "no_intervals" ? meta.sample + "_" + intervals.baseName : meta.sample
            intervals = intervals.baseName != "no_intervals" ? intervals : []
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
        ch_versions = ch_versions.mix(FREEBAYES.out.versions)


        if(no_intervals){
            TABIX_FREEBAYES(FREEBAYES.out.vcf)
            freebayes_vcf_gz_tbi = FREEBAYES.out.vcf.join(TABIX_FREEBAYES.out.tbi)
            ch_versions = ch_versions.mix(TABIX_FREEBAYES.out.versions)
        }else{
            BGZIP_FREEBAYES(FREEBAYES.out.vcf)
            freebayes_vcf_to_concat = BGZIP_FREEBAYES.out.vcf.groupTuple(size: num_intervals)

            CONCAT_VCF_FREEBAYES(freebayes_vcf_to_concat,fasta_fai, intervals_bed_combine_gz)
            freebayes_vcf_gz_tbi = CONCAT_VCF_FREEBAYES.out.vcf

            ch_versions = ch_versions.mix(BGZIP_FREEBAYES.out.versions)
            ch_versions = ch_versions.mix(CONCAT_VCF_FREEBAYES.out.versions)
        }
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

        ch_versions = ch_versions.mix(HAPLOTYPECALLER.out.versions)

        if(no_intervals){
            TABIX_HAPLOTYPECALLER(HAPLOTYPECALLER.out.vcf)
            haplotypecaller_gvcf_gz_tbi = HAPLOTYPECALLER.out.vcf.join(TABIX_HAPLOTYPECALLER.out.tbi)
            ch_versions = ch_versions.mix(TABIX_HAPLOTYPECALLER.out.versions)
        }else{
            BGZIP_HAPLOTYPECALLER(HAPLOTYPECALLER.out.vcf)

            haplotypecaller_gvcf_to_concat = BGZIP_HAPLOTYPECALLER.out.vcf.groupTuple(size: num_intervals)
            CONCAT_VCF_HAPLOTYPECALLER(haplotypecaller_gvcf_to_concat, fasta_fai, intervals_bed_combine_gz)
            haplotypecaller_gvcf_gz_tbi = CONCAT_VCF_HAPLOTYPECALLER.out.vcf

            ch_versions = ch_versions.mix(BGZIP_HAPLOTYPECALLER.out.versions)
            ch_versions = ch_versions.mix(CONCAT_VCF_HAPLOTYPECALLER.out.versions)
        }

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
        //TODO: Research if splitting by intervals is ok, we pretend for now it is fine. Seems to be the consensus on upstream modules implementaiton too

        MANTA_GERMLINE(
            cram_recalibrated_intervals_gz_tbi,
            fasta,
            fasta_fai
        )

        ch_versions = ch_versions.mix(MANTA_GERMLINE.out.versions)

        if(no_intervals){
            manta_candidate_small_indels_vcf_tbi = MANTA_GERMLINE.out.candidate_small_indels_vcf.join(MANTA_GERMLINE.out.candidate_small_indels_vcf_tbi)
            manta_candidate_sv_vcf_tbi           = MANTA_GERMLINE.out.candidate_sv_vcf.join(MANTA_GERMLINE.out.candidate_sv_vcf_tbi)
            manta_diploid_sv_vcf_tbi             = MANTA_GERMLINE.out.diploid_sv_vcf.join(MANTA_GERMLINE.out.diploid_sv_vcf)
        }else{

            BGZIP_MANTA_SV(MANTA_GERMLINE.out.candidate_small_indels_vcf)
            BGZIP_MANTA_SMALL_INDELS(MANTA_GERMLINE.out.candidate_sv_vcf)
            BGZIP_MANTA_DIPLOID(MANTA_GERMLINE.out.diploid_sv_vcf)

            manta_sv_vcf_to_concat = BGZIP_MANTA_SV.out.vcf.groupTuple(size: num_intervals)
            manta_small_indels_vcf_to_concat = BGZIP_MANTA_SMALL_INDELS.out.vcf.groupTuple(size: num_intervals)
            manta_diploid_vcf_to_concat = BGZIP_MANTA_DIPLOID.out.vcf.groupTuple(size: num_intervals)

            CONCAT_VCF_MANTA_SV(manta_sv_vcf_to_concat, fasta_fai, intervals_bed_combine_gz)
            CONCAT_VCF_MANTA_SMALL_INDELS(manta_small_indels_vcf_to_concat,fasta_fai, intervals_bed_combine_gz)
            CONCAT_VCF_MANTA_DIPLOID(manta_diploid_vcf_to_concat, fasta_fai, intervals_bed_combine_gz)

            manta_candidate_small_indels_vcf_tbi = CONCAT_VCF_MANTA_SV.out.vcf
            manta_candidate_sv_vcf_tbi           = CONCAT_VCF_MANTA_SMALL_INDELS.out.vcf
            manta_diploid_sv_vcf_tbi             = CONCAT_VCF_MANTA_DIPLOID.out.vcf

            ch_versions = ch_versions.mix(BGZIP_MANTA_SV.out.versions)
            ch_versions = ch_versions.mix(BGZIP_MANTA_SMALL_INDELS.out.versions)
            ch_versions = ch_versions.mix(BGZIP_MANTA_DIPLOID.out.versions)

            ch_versions = ch_versions.mix(CONCAT_VCF_MANTA_SV.out.versions)
            ch_versions = ch_versions.mix(CONCAT_VCF_MANTA_SMALL_INDELS.out.versions)
            ch_versions = ch_versions.mix(CONCAT_VCF_MANTA_DIPLOID.out.versions)

        }
    }

    if (tools.contains('strelka')) {
        //TODO: Research if splitting by intervals is ok, no reply on issue,  we pretend for now it is fine. Seems to be the consensus on upstream modules implementaiton too

        STRELKA_GERMLINE(
            cram_recalibrated_intervals_gz_tbi,
            fasta,
            fasta_fai
            )

        ch_versions = ch_versions.mix(STRELKA_GERMLINE.out.versions)

        if(no_intervals){
            strelka_vcf_gz_tbi = STRELKA_GERMLINE.out.vcf.join(STRELKA_GERMLINE.out.vcf_tbi)
        }else{
            BGZIP_STRELKA(STRELKA_GERMLINE.out.vcf)
            strelka_vcf_to_concat = BGZIP_STRELKA.out.vcf.groupTuple(size: num_intervals)

            CONCAT_VCF_STRELKA(strelka_vcf_to_concat,fasta_fai, intervals_bed_combine_gz)
            strelka_vcf_gz_tbi = CONCAT_VCF_STRELKA.out.vcf

            ch_versions = ch_versions.mix(BGZIP_STRELKA.out.versions)
            ch_versions = ch_versions.mix(CONCAT_VCF_STRELKA.out.versions)
        }
    }

    if (tools.contains('tiddit')){
        //TODO: Update tiddit on bioconda, the current version does not support cram usage, needs newest version:
        // https://github.com/SciLifeLab/TIDDIT/issues/82#issuecomment-1022103264
        // Issue opened, either this week or end of february

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
    manta_candidate_small_indels_vcf_tbi
    manta_candidate_sv_vcf_tbi
    manta_diploid_sv_vcf_tbi
    strelka_vcf_gz_tbi

    versions = ch_versions
}
