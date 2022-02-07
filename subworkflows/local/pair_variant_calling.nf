//
// PAIRED VARIANT CALLING
//
include { BGZIP as BGZIP_MANTA_SMALL_INDELS           } from '../../modules/local/bgzip'
include { BGZIP as BGZIP_MANTA_SV                     } from '../../modules/local/bgzip'
include { BGZIP as BGZIP_MANTA_DIPLOID                } from '../../modules/local/bgzip'
include { BGZIP as BGZIP_MANTA_SOMATIC                  } from '../../modules/local/bgzip'
include { BGZIP as BGZIP_STRELKA                      } from '../../modules/local/bgzip'
include { BGZIP as BGZIP_STRELKA_BP                      } from '../../modules/local/bgzip'
include { CONCAT_VCF as CONCAT_VCF_MANTA_SMALL_INDELS } from '../../modules/local/concat_vcf/main'
include { CONCAT_VCF as CONCAT_VCF_MANTA_SV           } from '../../modules/local/concat_vcf/main'
include { CONCAT_VCF as CONCAT_VCF_MANTA_DIPLOID      } from '../../modules/local/concat_vcf/main'
include { CONCAT_VCF as CONCAT_VCF_MANTA_SOMATIC        } from '../../modules/local/concat_vcf/main'
include { CONCAT_VCF as CONCAT_VCF_STRELKA            } from '../../modules/local/concat_vcf/main'
include { CONCAT_VCF as CONCAT_VCF_STRELKA_BP            } from '../../modules/local/concat_vcf/main'
include { GATK_TUMOR_NORMAL_SOMATIC_VARIANT_CALLING     } from '../../subworkflows/nf-core/gatk_tumor_normal_somatic_variant_calling/main'
include { MANTA_SOMATIC                                 } from '../../modules/nf-core/modules/manta/somatic/main'
include { MSISENSORPRO_MSI_SOMATIC                      } from '../../modules/nf-core/modules/msisensorpro/msi_somatic/main'
include { STRELKA_SOMATIC                               } from '../../modules/nf-core/modules/strelka/somatic/main'
include { STRELKA_SOMATIC as STRELKA_BP                 } from '../../modules/nf-core/modules/strelka/somatic/main'

workflow PAIR_VARIANT_CALLING {
    take:
        tools
        cram_pair             // channel: [mandatory] cram
        dbsnp                 // channel: [mandatory] dbsnp
        dbsnp_tbi             // channel: [mandatory] dbsnp_tbi
        dict                  // channel: [mandatory] dict
        fasta                 // channel: [mandatory] fasta
        fasta_fai             // channel: [mandatory] fasta_fai
        intervals               // channel: [mandatory] intervals/target regions
        intervals_bed_gz_tbi    // channel: [mandatory] intervals/target regions index zipped and indexed
        intervals_bed_combined_gz_tbi    // channel: [mandatory] intervals/target regions index zipped and indexed
        intervals_bed_combine_gz            // channel: [mandatory] intervals/target regions index zipped and indexed in one file
        num_intervals                       // val: number of intervals that are used to parallelize exection, either based on capture kit or GATK recommended for WGS
        no_intervals
        msisensorpro_scan     // channel: [optional]  msisensorpro_scan
        germline_resource     // channel: [optional]  germline_resource
        germline_resource_tbi // channel: [optional]  germline_resource_tbi
        panel_of_normals      // channel: [optional]  panel_of_normals
        panel_of_normals_tbi  // channel: [optional]  panel_of_normals_tbi

    main:

    if(!tools) tools = ""

    ch_versions          = Channel.empty()
    manta_vcf            = Channel.empty()
    strelka_vcf          = Channel.empty()

    cram_pair.combine(intervals)
        .map{ meta, normal_cram, normal_crai, tumor_cram, tumor_crai, intervals ->
            new_meta = meta.clone()
            new_meta.id = intervals.baseName != "no_intervals" ? meta.tumor_id + "_vs_" + meta.normal_id + "_" + intervals.baseName : meta.sample
            intervals = intervals.baseName != "no_intervals" ? intervals : []
            [new_meta, normal_cram, normal_crai, tumor_cram, tumor_crai, intervals]
        }.set{cram_pair_intervals}

    cram_pair.combine(intervals_bed_gz_tbi)
        .map{ meta, normal_cram, normal_crai, tumor_cram, tumor_crai, bed, tbi ->
            new_meta = meta.clone()
            new_meta.id = bed.simpleName != "no_intervals" ? meta.tumor_id + "_vs_" + meta.normal_id + "_" + bed.simpleName : meta.sample
            bed = bed.simpleName != "no_intervals" ? bed : []
            tbi = tbi.simpleName != "no_intervals" ? tbi : []

            [new_meta, normal_cram, normal_crai, tumor_cram, tumor_crai, bed, tbi]
        }.set{cram_pair_intervals_gz_tbi}

    if (tools.contains('manta')) {
        MANTA_SOMATIC(
            cram_pair_intervals_gz_tbi,
            fasta,
            fasta_fai,
        )

        ch_versions = ch_versions.mix(MANTA_SOMATIC.out.versions)

        if(no_intervals){
            manta_candidate_small_indels_vcf_tbi = MANTA_SOMATIC.out.candidate_small_indels_vcf.join(MANTA_SOMATIC.out.candidate_small_indels_vcf_tbi)
            manta_candidate_sv_vcf_tbi           = MANTA_SOMATIC.out.candidate_sv_vcf.join(MANTA_SOMATIC.out.candidate_sv_vcf_tbi)
            manta_diploid_sv_vcf_tbi             = MANTA_SOMATIC.out.diploid_sv_vcf.join(MANTA_SOMATIC.out.diploid_sv_vcf)
            manta_somatic_sv_vcf_tbi             = MANTA_SOMATIC.out.somatic_sv_vcf.join(MANTA_SOMATIC.out.somatic_sv_vcf)
        }else{

            BGZIP_MANTA_SV(MANTA_SOMATIC.out.candidate_small_indels_vcf)
            BGZIP_MANTA_SMALL_INDELS(MANTA_SOMATIC.out.candidate_sv_vcf)
            BGZIP_MANTA_DIPLOID(MANTA_SOMATIC.out.diploid_sv_vcf)
            BGZIP_MANTA_SOMATIC(MANTA_SOMATIC.out.somatic_sv_vcf)

            BGZIP_MANTA_SV.out.vcf.map{ meta, vcf ->
                new_meta = meta.clone()
                new_meta.id = new_meta.tumor_id + "_vs_" + new_meta.normal_id
                [new_meta, vcf]
            }.groupTuple(size: num_intervals)
            .set{manta_sv_vcf_to_concat}

            BGZIP_MANTA_SMALL_INDELS.out.vcf.map{ meta, vcf ->
                new_meta = meta.clone()
                new_meta.id = new_meta.tumor_id + "_vs_" + new_meta.normal_id
                [new_meta, vcf]
            }.groupTuple(size: num_intervals)
            .set{manta_small_indels_vcf_to_concat}

            BGZIP_MANTA_DIPLOID.out.vcf.map{ meta, vcf ->
                new_meta = meta.clone()
                new_meta.id = new_meta.tumor_id + "_vs_" + new_meta.normal_id
                [new_meta, vcf]
            }.groupTuple(size: num_intervals)
            .set{manta_diploid_vcf_to_concat}

            BGZIP_MANTA_SOMATIC.out.vcf.map{ meta, vcf ->
                new_meta = meta.clone()
                new_meta.id = new_meta.tumor_id + "_vs_" + new_meta.normal_id
                [new_meta, vcf]
            }.groupTuple(size: num_intervals)
            .set{manta_somatic_sv_vcf_to_concat}

            CONCAT_VCF_MANTA_SV(manta_sv_vcf_to_concat, fasta_fai, intervals_bed_combine_gz)
            CONCAT_VCF_MANTA_SMALL_INDELS(manta_small_indels_vcf_to_concat,fasta_fai, intervals_bed_combine_gz)
            CONCAT_VCF_MANTA_DIPLOID(manta_diploid_vcf_to_concat, fasta_fai, intervals_bed_combine_gz)
            CONCAT_VCF_MANTA_SOMATIC(manta_somatic_sv_vcf_to_concat, fasta_fai, intervals_bed_combine_gz)

            manta_candidate_small_indels_vcf_tbi = CONCAT_VCF_MANTA_SV.out.vcf
            manta_candidate_sv_vcf_tbi           = CONCAT_VCF_MANTA_SMALL_INDELS.out.vcf
            manta_diploid_sv_vcf_tbi             = CONCAT_VCF_MANTA_DIPLOID.out.vcf
            manta_somatic_sv_vcf_tbi             = CONCAT_VCF_MANTA_SOMATIC.out.vcf

            ch_versions = ch_versions.mix(BGZIP_MANTA_SV.out.versions)
            ch_versions = ch_versions.mix(BGZIP_MANTA_SMALL_INDELS.out.versions)
            ch_versions = ch_versions.mix(BGZIP_MANTA_DIPLOID.out.versions)
            ch_versions = ch_versions.mix(BGZIP_MANTA_SOMATIC.out.versions)

            ch_versions = ch_versions.mix(CONCAT_VCF_MANTA_SV.out.versions)
            ch_versions = ch_versions.mix(CONCAT_VCF_MANTA_SMALL_INDELS.out.versions)
            ch_versions = ch_versions.mix(CONCAT_VCF_MANTA_DIPLOID.out.versions)
            ch_versions = ch_versions.mix(CONCAT_VCF_MANTA_SOMATIC.out.versions)

        }

        //manta_vcf = manta_candidate_small_indels_vcf.mix(manta_candidate_sv_vcf,manta_diploid_sv_vcf,manta_somatic_sv_vcf)

        if (tools.contains('strelka')) {
            cram_pair.join(manta_somatic_sv_vcf_tbi).combine(intervals_bed_gz_tbi)
            .map{ meta, normal_cram, normal_crai, tumor_cram, tumor_crai, manta_vcf, manta_tbi, bed, tbi ->
                new_meta = meta.clone()
                new_meta.id = bed.simpleName != "no_intervals" ? meta.tumor_id + "_vs_" + meta.normal_id + "_" + bed.simpleName : meta.sample
                bed = bed.simpleName != "no_intervals" ? bed : []
                tbi = tbi.simpleName != "no_intervals" ? tbi : []

                [new_meta, normal_cram, normal_crai, tumor_cram, tumor_crai, manta_vcf, manta_tbi, bed, tbi]
            }.set{cram_pair_strelka}

            STRELKA_BP(
                cram_pair_strelka,
                fasta,
                fasta_fai
                )

            if(no_intervals){
                strelka_snvs_vcf_gz_tbi = STRELKA_BP.out.vcf_snvs.join(STRELKA_BP.out.vcf_snvs_tbi)
            }else{
                BGZIP_STRELKA_BP(STRELKA_BP.out.vcf_snvs)

                BGZIP_STRELKA_BP.out.vcf.map{ meta, vcf ->
                    new_meta = meta.clone()
                    new_meta.id = new_meta.tumor_id + "_vs_" + new_meta.normal_id
                    [new_meta, vcf]
                }.groupTuple(size: num_intervals)
                .set{strelka_vcf_to_concat}

                CONCAT_VCF_STRELKA_BP(strelka_vcf_to_concat,fasta_fai, intervals_bed_combine_gz)
                strelka_vcf_gz_tbi = CONCAT_VCF_STRELKA_BP.out.vcf

                ch_versions = ch_versions.mix(BGZIP_STRELKA_BP.out.versions)
                ch_versions = ch_versions.mix(CONCAT_VCF_STRELKA_BP.out.versions)
            }
            //strelka_vcf = strelka_vcf.mix(strelka_indels_vcf,strelka_snvs_vcf)
        }
    }else{
        if (tools.contains('strelka')) {
            //TODO: research if multiple targets can be provided: waiting for reply
            cram_pair.combine(intervals_bed_gz_tbi)
            .map{ meta, normal_cram, normal_crai, tumor_cram, tumor_crai, bed, tbi ->
                new_meta = meta.clone()
                new_meta.id = bed.simpleName != "no_intervals" ? meta.tumor_id + "_vs_" + meta.normal_id + "_" + bed.simpleName : meta.sample
                bed = bed.simpleName != "no_intervals" ? bed : []
                tbi = tbi.simpleName != "no_intervals" ? tbi : []

                [new_meta, normal_cram, normal_crai, tumor_cram, tumor_crai, [], [], bed, tbi]
            }.set{cram_pair_strelka}

            STRELKA_SOMATIC(
                cram_pair_strelka,
                fasta,
                fasta_fai)

            ch_versions = ch_versions.mix(STRELKA_SOMATIC.out.versions)

            if(no_intervals){
                strelka_snvs_vcf_gz_tbi = STRELKA_SOMATIC.out.vcf_snvs.join(STRELKA_SOMATIC.out.vcf_snvs_tbi)
            }else{
                BGZIP_STRELKA(STRELKA_SOMATIC.out.vcf_snvs)

                BGZIP_STRELKA.out.vcf.map{ meta, vcf ->
                    new_meta = meta.clone()
                    new_meta.id = new_meta.tumor_id + "_vs_" + new_meta.normal_id
                    [new_meta, vcf]
                }.groupTuple(size: num_intervals)
                .set{strelka_vcf_to_concat}

                CONCAT_VCF_STRELKA(strelka_vcf_to_concat,fasta_fai, intervals_bed_combine_gz)
                strelka_vcf_gz_tbi = CONCAT_VCF_STRELKA.out.vcf

                ch_versions = ch_versions.mix(BGZIP_STRELKA.out.versions)
                ch_versions = ch_versions.mix(CONCAT_VCF_STRELKA.out.versions)
            }
        }
    }

    if (tools.contains('msisensorpro')) {
        cram_pair.view()
        fasta.view()
        msisensorpro_scan.view()
        MSISENSORPRO_MSI_SOMATIC(
            cram_pair,
            fasta,
            msisensorpro_scan)
        ch_versions = ch_versions.mix(MSISENSORPRO_MSI_SOMATIC.out.versions)
    }

    if (tools.contains('mutect2')){
        //need some feedback here from Gavin
        which_norm = []
        cram_pair_intervals.map{ meta, normal_cram, normal_crai, tumor_cram, tumor_crai, intervals -> [meta, normal_cram, normal_crai, tumor_cram, tumor_crai, intervals, which_norm]}.set{cram_pair_mutect2}

        // GATK_TUMOR_NORMAL_SOMATIC_VARIANT_CALLING(
        //     cram_pair_mutect2,
        //     fasta,
        //     fasta_fai,
        //     dict,
        //     germline_resource,
        //     germline_resource_tbi,
        //     panel_of_normals,
        //     panel_of_normals_tbi,
        //     no_intervals
        // )
        //ch_versions = ch_versions.mix(GATK_TUMOR_NORMAL_SOMATIC_VARIANT_CALLING.out.versions)
    }

    // if (tools.contains('tiddit')){
    // }

    emit:
    versions    = ch_versions
        // manta_vcf            = manta_vcf
        // strelka_vcf          = strelka_vcf
}
