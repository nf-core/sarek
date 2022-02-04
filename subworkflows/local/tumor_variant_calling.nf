//
// TUMOR VARIANT CALLING
// Should be only run on patients without normal sample
//


include { BGZIP as BGZIP_FREEBAYES                    } from '../../modules/local/bgzip'
include { BGZIP as BGZIP_MANTA_SMALL_INDELS           } from '../../modules/local/bgzip'
include { BGZIP as BGZIP_MANTA_SV                     } from '../../modules/local/bgzip'
include { BGZIP as BGZIP_MANTA_TUMOR                  } from '../../modules/local/bgzip'
include { BGZIP as BGZIP_STRELKA                      } from '../../modules/local/bgzip'
include { CONCAT_VCF as CONCAT_VCF_FREEBAYES          } from '../../modules/local/concat_vcf/main'
include { CONCAT_VCF as CONCAT_VCF_MANTA_SMALL_INDELS } from '../../modules/local/concat_vcf/main'
include { CONCAT_VCF as CONCAT_VCF_MANTA_SV           } from '../../modules/local/concat_vcf/main'
include { CONCAT_VCF as CONCAT_VCF_MANTA_TUMOR        } from '../../modules/local/concat_vcf/main'
include { CONCAT_VCF as CONCAT_VCF_STRELKA            } from '../../modules/local/concat_vcf/main'
include { FREEBAYES                                   } from '../../modules/nf-core/modules/freebayes/main'
include { GATK_TUMOR_ONLY_SOMATIC_VARIANT_CALLING     } from '../../subworkflows/nf-core/gatk_tumor_only_somatic_variant_calling/main'
include { MANTA_TUMORONLY                             } from '../../modules/nf-core/modules/manta/tumoronly/main'
include { STRELKA_GERMLINE as STRELKA_TUMORONLY       } from '../../modules/nf-core/modules/strelka/germline/main'

workflow TUMOR_ONLY_VARIANT_CALLING {
    take:
        tools                           // Mandatory, list of tools to apply
        cram_recalibrated               // channel: [mandatory] cram
        dbsnp                           // channel: [mandatory] dbsnp
        dbsnp_tbi                       // channel: [mandatory] dbsnp_tbi
        dict                            // channel: [mandatory] dict
        fasta                           // channel: [mandatory] fasta
        fasta_fai                       // channel: [mandatory] fasta_fai
        intervals                       // channel: [mandatory] intervals/target regions
        intervals_bed_gz_tbi            // channel: [mandatory] intervals/target regions index zipped and indexed
        intervals_bed_combine_gz_tbi    // channel: [mandatory] intervals/target regions index zipped and indexed
        intervals_bed_combine_gz        // channel: [mandatory] intervals/target regions index zipped and indexed in one file
        num_intervals                   // val: number of intervals that are used to parallelize exection, either based on capture kit or GATK recommended for WGS
        no_intervals
        germline_resource
        germline_resource_tbi   // channel
        panel_of_normals
        panel_of_normals_tbi


    main:

    ch_versions                             = Channel.empty()
    freebayes_vcf_gz_tbi                    = Channel.empty()
    manta_candidate_small_indels_vcf_tbi    = Channel.empty()
    manta_candidate_sv_vcf_tbi              = Channel.empty()
    manta_tumor_sv_vcf_tbi                  = Channel.empty()
    mutect2_vcf_gz_tbi                      = Channel.empty()
    strelka_vcf_gz_tbi                      = Channel.empty()

    cram_recalibrated.combine(intervals)
        .map{ meta, cram, crai, intervals ->
            new_meta = meta.clone()
            new_meta.id = intervals.baseName != "no_intervals" ? meta.sample + "_" + intervals.baseName : meta.sample
            intervals = intervals.baseName != "no_intervals" ? intervals : []
            [new_meta, cram, crai, intervals]
        }.set{cram_recalibrated_intervals}


    cram_recalibrated.combine(intervals_bed_gz_tbi)
        .map{ meta, cram, crai, bed, tbi ->
            new_meta = meta.clone()
            new_meta.id = bed.simpleName != "no_intervals" ? meta.sample + "_" + bed.simpleName : meta.sample
            bed = bed.simpleName != "no_intervals" ? bed : []
            tbi = tbi.simpleName != "no_intervals" ? tbi : []

            [new_meta, cram, crai, bed, tbi]
        }.set{cram_recalibrated_intervals_gz_tbi}

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
            BGZIP_FREEBAYES.out.vcf.map{ meta, vcf ->
                new_meta = meta.clone()
                new_meta.id = new_meta.sample
                [new_meta, vcf]
            }.groupTuple(size: num_intervals)
            .set{freebayes_vcf_to_concat}

            CONCAT_VCF_FREEBAYES(freebayes_vcf_to_concat,fasta_fai, intervals_bed_combine_gz)
            freebayes_vcf_gz_tbi = CONCAT_VCF_FREEBAYES.out.vcf

            ch_versions = ch_versions.mix(BGZIP_FREEBAYES.out.versions)
            ch_versions = ch_versions.mix(CONCAT_VCF_FREEBAYES.out.versions)
        }
    }

    if (tools.contains('mutect2')) {

        which_norm = []
        cram_recalibrated_intervals.map{ meta, cram, crai, intervals -> [meta, cram, crai, intervals, which_norm]}.set{cram_recalibrated_mutect2}
        GATK_TUMOR_ONLY_SOMATIC_VARIANT_CALLING(
            cram_recalibrated_mutect2,
            fasta,
            fasta_fai,
            dict,
            germline_resource,
            germline_resource_tbi,
            panel_of_normals,
            panel_of_normals_tbi,
            num_intervals,
            no_intervals,
            intervals_bed_combine_gz
        )

        ch_versions = ch_versions.mix(GATK_TUMOR_ONLY_SOMATIC_VARIANT_CALLING.out.versions)
    }

    if (tools.contains('manta')){
        //TODO: Research if splitting by intervals is ok, we pretend for now it is fine. Seems to be the consensus on upstream modules implementaiton too

        MANTA_TUMORONLY(
            cram_recalibrated_intervals_gz_tbi,
            fasta,
            fasta_fai
        )

        ch_versions = ch_versions.mix(MANTA_TUMORONLY.out.versions)

        if(no_intervals){
            manta_candidate_small_indels_vcf_tbi = MANTA_TUMORONLY.out.candidate_small_indels_vcf.join(MANTA_TUMORONLY.out.candidate_small_indels_vcf_tbi)
            manta_candidate_sv_vcf_tbi           = MANTA_TUMORONLY.out.candidate_sv_vcf.join(MANTA_TUMORONLY.out.candidate_sv_vcf_tbi)
            manta_tumor_sv_vcf_tbi               = MANTA_TUMORONLY.out.tumor_sv_vcf.join(MANTA_TUMORONLY.out.tumor_sv_vcf)
        }else{

            BGZIP_MANTA_SV(MANTA_TUMORONLY.out.candidate_small_indels_vcf)
            BGZIP_MANTA_SMALL_INDELS(MANTA_TUMORONLY.out.candidate_sv_vcf)
            BGZIP_MANTA_TUMOR(MANTA_TUMORONLY.out.tumor_sv_vcf)

            BGZIP_MANTA_SV.out.vcf.map{ meta, vcf ->
                new_meta = meta.clone()
                new_meta.id = new_meta.sample
                [new_meta, vcf]
            }.groupTuple(size: num_intervals)
            .set{manta_sv_vcf_to_concat}

            BGZIP_MANTA_SMALL_INDELS.out.vcf.map{ meta, vcf ->
                new_meta = meta.clone()
                new_meta.id = new_meta.sample
                [new_meta, vcf]
            }.groupTuple(size: num_intervals)
            .set{manta_small_indels_vcf_to_concat}

            BGZIP_MANTA_TUMOR.out.vcf.map{ meta, vcf ->
                new_meta = meta.clone()
                new_meta.id = new_meta.sample
                [new_meta, vcf]
            }.groupTuple(size: num_intervals)
            .set{manta_tumor_sv_vcf_to_concat}

            CONCAT_VCF_MANTA_SV(manta_sv_vcf_to_concat, fasta_fai, intervals_bed_combine_gz)
            CONCAT_VCF_MANTA_SMALL_INDELS(manta_small_indels_vcf_to_concat,fasta_fai, intervals_bed_combine_gz)
            CONCAT_VCF_MANTA_TUMOR(manta_tumor_sv_vcf_to_concat, fasta_fai, intervals_bed_combine_gz)

            manta_candidate_small_indels_vcf_tbi = CONCAT_VCF_MANTA_SV.out.vcf
            manta_candidate_sv_vcf_tbi           = CONCAT_VCF_MANTA_SMALL_INDELS.out.vcf
            manta_tumor_sv_vcf_tbi               = CONCAT_VCF_MANTA_TUMOR.out.vcf

            ch_versions = ch_versions.mix(BGZIP_MANTA_SV.out.versions)
            ch_versions = ch_versions.mix(BGZIP_MANTA_SMALL_INDELS.out.versions)
            ch_versions = ch_versions.mix(BGZIP_MANTA_TUMOR.out.versions)

            ch_versions = ch_versions.mix(CONCAT_VCF_MANTA_SV.out.versions)
            ch_versions = ch_versions.mix(CONCAT_VCF_MANTA_SMALL_INDELS.out.versions)
            ch_versions = ch_versions.mix(CONCAT_VCF_MANTA_TUMOR.out.versions)

        }
    }

    if (tools.contains('strelka')) {
        //TODO: research if multiple targets can be provided: waiting for reply

        STRELKA_TUMORONLY(
            cram_recalibrated_intervals_gz_tbi,
            fasta,
            fasta_fai
            )

        ch_versions = ch_versions.mix(STRELKA_TUMORONLY.out.versions)

        if(no_intervals){
            strelka_vcf_gz_tbi = STRELKA_TUMORONLY.out.vcf.join(STRELKA_TUMORONLY.out.vcf_tbi)
        }else{
            BGZIP_STRELKA(STRELKA_TUMORONLY.out.vcf)

            BGZIP_STRELKA.out.vcf.map{ meta, vcf ->
                new_meta = meta.clone()
                new_meta.id = new_meta.sample
                [new_meta, vcf]
            }.groupTuple(size: num_intervals)
            .set{strelka_vcf_to_concat}

            CONCAT_VCF_STRELKA(strelka_vcf_to_concat,fasta_fai, intervals_bed_combine_gz)
            strelka_vcf_gz_tbi = CONCAT_VCF_STRELKA.out.vcf

            ch_versions = ch_versions.mix(BGZIP_STRELKA.out.versions)
            ch_versions = ch_versions.mix(CONCAT_VCF_STRELKA.out.versions)
        }
    }


    if (tools.contains('tiddit')){
    }

    emit:
    ch_versions
    freebayes_vcf_gz_tbi
    manta_candidate_small_indels_vcf_tbi
    manta_candidate_sv_vcf_tbi
    manta_tumor_sv_vcf_tbi
    mutect2_vcf_gz_tbi
    strelka_vcf_gz_tbi
}
