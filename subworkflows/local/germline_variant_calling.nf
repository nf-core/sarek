//
// GERMLINE VARIANT CALLING
//

include { GATK_JOINT_GERMLINE_VARIANT_CALLING        } from '../../subworkflows/nf-core/joint_germline_variant_calling/main'
include { DEEPVARIANT                                } from '../../modules/nf-core/modules/deepvariant/main'
include { FREEBAYES                                  } from '../../modules/nf-core/modules/freebayes/main'
include { GATK4_HAPLOTYPECALLER as HAPLOTYPECALLER   } from '../../modules/nf-core/modules/gatk4/haplotypecaller/main'
include { MANTA_GERMLINE                             } from '../../modules/nf-core/modules/manta/germline/main'
include { STRELKA_GERMLINE                           } from '../../modules/nf-core/modules/strelka/germline/main'
include { TIDDIT_SV                                  } from '../../modules/nf-core/modules/tiddit/sv/main'
include { TABIX_BGZIPTABIX as TABIX_BGZIP_TIDDIT_SV  } from '../../modules/nf-core/modules/tabix/bgziptabix/main'

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
        //target_bed        // channel: [optional]  target_bed
        //target_bed_gz_tbi // channel: [optional]  target_bed_gz_tbi

    main:

    ch_versions = Channel.empty()

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


    haplotypecaller_gvcf = Channel.empty()
    haplotypecaller_vcf  = Channel.empty()
    strelka_vcf          = Channel.empty()

    if (tools.contains('haplotypecaller')) {

        //TODO: merge parallelized vcfs, index them all
        //TODO: Pass over dbsnp/knwon_indels?
        HAPLOTYPECALLER(
                cram_recalibrated_intervals,
                fasta,
                fasta_fai,
                dict,
                dbsnp,
                dbsnp_tbi
        )
        ch_versions = ch_versions.mix(HAPLOTYPECALLER.out.versions)

        haplotypecaller_vcf_gz_tbi = HAPLOTYPECALLER.out.vcf.join(HAPLOTYPECALLER.out.tbi)

        if(joint_germline){
            run_haplotypecaller = false
            run_vqsr            = true
            //Waiting on some feedback from gavin
            // GATK_JOINT_GERMLINE_VARIANT_CALLING(
            //     haplotypecaller_vcf_gz_tbi,
            //     run_haplotypecaller,
            //     run_vqsr,
            //     fasta,
            //     fasta_fai,
            //     dict,
            //     //TODO: replace with known_sites?
            //     dbsnp,
            //     dbsnp_tbi,
            //     "joined",

            // )
        }

    }

    if (tools.contains('deepvariant')) {
        //TODO: merge parallelized vcfs, index them all
        //TODO: research if multiple targets can be provided
        //TODO: Pass over dbsnp/knwon_indels?
        DEEPVARIANT(
            cram_recalibrated_intervals,
            fasta,
            fasta_fai)

        deepvariant_vcf = DEEPVARIANT.out.vcf
        deepvariant_gvcf = DEEPVARIANT.out.gvcf
    }

    if (tools.contains('freebayes')){
        //TODO: merge parallelized vcfs, index them all
        //TODO: research if multiple targets can be provided
        //TODO: Pass over dbsnp/knwon_indels?
        cram_recalibrated.combine(intervals)
        .map{ meta, cram, crai, intervals ->
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
        freebayes_vcf_gz = FREEBAYES.out.vcf
        ch_versions = ch_versions.mix(FREEBAYES.out.versions)
    }

    if (tools.contains('manta')){
        //TODO: test data not running
        //TODO: merge parallelized vcfs, index them all
        //TODO: Pass over dbsnp/knwon_indels?
        cram_recalibrated
        .map{ meta, cram, crai ->
            new_meta = meta.clone()
            //new_meta.id = meta.sample + "_" + intervals.simpleName
            [new_meta, cram, crai, [], []]
        }.set{cram_recalibrated_manta}

        MANTA_GERMLINE(
            cram_recalibrated_manta,
            fasta,
            fasta_fai
        )
        manta_candidate_small_indels_vcf_tbi = MANTA_GERMLINE.out.candidate_small_indels_vcf.join(MANTA_GERMLINE.out.candidate_small_indels_vcf_tbi)
        manta_candidate_sv_vcf_tbi = MANTA_GERMLINE.out.candidate_sv_vcf.join(MANTA_GERMLINE.out.candidate_sv_vcf_tbi)
        manta_diploid_sv_vcf_tbi = MANTA_GERMLINE.out.diploid_sv_vcf.join(MANTA_GERMLINE.out.diploid_sv_vcf)

        ch_versions = ch_versions.mix(MANTA_GERMLINE.out.versions)
    }

    if (tools.contains('strelka')) {
        //TODO: research if multiple targets can be provided
        //TODO: merge parallelized vcfs, index them all
        //TODO: Pass over dbsnp/knwon_indels?

        STRELKA_GERMLINE(
            cram_recalibrated_intervals_gz_tbi,
            fasta,
            fasta_fai
            )

        strelka_vcf_tbi = STRELKA_GERMLINE.out.vcf.join(STRELKA_GERMLINE.out.vcf_tbi)
        strelka_genome_vcf_tbi = STRELKA_GERMLINE.out.genome_vcf.join(STRELKA_GERMLINE.out.genome_vcf_tbi)

        ch_versions = ch_versions.mix(STRELKA_GERMLINE.out.versions)

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

        ch_versions = ch_versions.mix(TIDDIT_SV.out.versions)
    }

    emit:
    strelka_vcf_tbi        = Channel.empty()
    strelka_genome_vcf_tbi = Channel.empty()
}
