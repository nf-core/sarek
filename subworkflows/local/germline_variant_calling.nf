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
include { BGZIP as BGZIP_STRELKA_GENOME               } from '../../modules/local/bgzip'
include { CONCAT_VCF as CONCAT_VCF_DEEPVARIANT        } from '../../modules/local/concat_vcf/main'
include { CONCAT_VCF as CONCAT_GVCF_DEEPVARIANT       } from '../../modules/local/concat_vcf/main'
include { CONCAT_VCF as CONCAT_VCF_FREEBAYES          } from '../../modules/local/concat_vcf/main'
include { CONCAT_VCF as CONCAT_VCF_HAPLOTYPECALLER    } from '../../modules/local/concat_vcf/main'
include { CONCAT_VCF as CONCAT_VCF_MANTA_SMALL_INDELS } from '../../modules/local/concat_vcf/main'
include { CONCAT_VCF as CONCAT_VCF_MANTA_SV           } from '../../modules/local/concat_vcf/main'
include { CONCAT_VCF as CONCAT_VCF_MANTA_DIPLOID      } from '../../modules/local/concat_vcf/main'
include { CONCAT_VCF as CONCAT_VCF_STRELKA            } from '../../modules/local/concat_vcf/main'
include { CONCAT_VCF as CONCAT_VCF_STRELKA_GENOME     } from '../../modules/local/concat_vcf/main'
include { DEEPVARIANT                                 } from '../../modules/nf-core/modules/deepvariant/main'
include { FREEBAYES                                   } from '../../modules/nf-core/modules/freebayes/main'
include { GATK_JOINT_GERMLINE_VARIANT_CALLING         } from '../../subworkflows/nf-core/joint_germline_variant_calling/main'
include { GATK4_HAPLOTYPECALLER as HAPLOTYPECALLER    } from '../../modules/nf-core/modules/gatk4/haplotypecaller/main'
include { GATK4_GENOTYPEGVCFS as GENOTYPEGVCFS        } from '../../modules/nf-core/modules/gatk4/genotypegvcfs/main'
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
        intervals_bed_combine_gz            // channel: [mandatory] intervals/target regions index zipped and indexed in one file
        num_intervals                       // val: number of intervals that are used to parallelize exection, either based on capture kit or GATK recommended for WGS
        no_intervals
        joint_germline                      // val: true/false on whether to run joint_germline calling, only works in combination with haplotypecaller at the moment

    main:

    if(!tools) tools = ""

    ch_versions = Channel.empty()

    deepvariant_vcf       = Channel.empty()
    freebayes_vcf         = Channel.empty()
    haplotypecaller_gvcf  = Channel.empty()
    genotypegvcfs_vcf     = Channel.empty()
    manta_vcf             = Channel.empty()
    strelka_vcf           = Channel.empty()

    cram_recalibrated.combine(intervals).map{ meta, cram, crai, intervals ->
            sample = meta.sample
            new_intervals = intervals.baseName != "no_intervals" ? intervals : []
            id = new_intervals ? sample + "_" + new_intervals.baseName : sample
            [[ id: id, sample: meta.sample, gender: meta.gender, status: meta.status, patient: meta.patient ], cram, crai, new_intervals]
        }.set{cram_recalibrated_intervals}

    cram_recalibrated.combine(intervals_bed_gz_tbi)
        .map{ meta, cram, crai, bed, tbi ->
            sample = meta.sample
            new_bed = bed.simpleName != "no_intervals" ? bed : []
            new_tbi = tbi.simpleName != "no_intervals" ? tbi : []
            id = new_bed ? sample + "_" + new_bed.simpleName : sample
            new_meta = [ id: id, sample: meta.sample, gender: meta.gender, status: meta.status, patient: meta.patient ]
            [new_meta, cram, crai, new_bed, new_tbi]
        }.set{cram_recalibrated_intervals_gz_tbi}

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

            deepvariant_vcf_gz = DEEPVARIANT.out.vcf
            deepvariant_gvcf_gz = DEEPVARIANT.out.gvcf

            ch_versions = ch_versions.mix(TABIX_DEEPVARIANT_VCF.out.versions)
            ch_versions = ch_versions.mix(TABIX_DEEPVARIANT_GVCF.out.versions)
        }else{
            BGZIP_DEEPVARIANT_VCF(DEEPVARIANT.out.vcf)
            BGZIP_DEEPVARIANT_GVCF(DEEPVARIANT.out.gvcf)

            BGZIP_DEEPVARIANT_VCF.out.vcf.map{ meta, vcf ->
                new_meta = meta.clone()
                new_meta.id = new_meta.sample
                [new_meta, vcf]
            }.groupTuple(size: num_intervals)
            .set{deepvariant_vcf_to_concat}

            BGZIP_DEEPVARIANT_GVCF.out.vcf.map{ meta, vcf ->
                new_meta = meta.clone()
                new_meta.id = new_meta.sample
                [new_meta, vcf]
            }.groupTuple(size: num_intervals)
            .set{deepvariant_gvcf_to_concat}

            CONCAT_VCF_DEEPVARIANT(deepvariant_vcf_to_concat, fasta_fai, intervals_bed_combine_gz)
            CONCAT_GVCF_DEEPVARIANT(deepvariant_gvcf_to_concat, fasta_fai, intervals_bed_combine_gz)

            deepvariant_vcf_gz = CONCAT_VCF_DEEPVARIANT.out.vcf
            deepvariant_gvcf_gz = CONCAT_GVCF_DEEPVARIANT.out.vcf

            ch_versions = ch_versions.mix(BGZIP_DEEPVARIANT_VCF.out.versions)
            ch_versions = ch_versions.mix(CONCAT_VCF_DEEPVARIANT.out.versions)
        }

        deepvariant_vcf = deepvariant_vcf.mix(deepvariant_vcf_gz, deepvariant_gvcf_gz)
    }

    if (tools.contains('freebayes')){

        cram_recalibrated.combine(intervals).map{ meta, cram, crai, intervals ->
            sample = meta.sample
            new_intervals = intervals.baseName != "no_intervals" ? intervals : []
            id = new_intervals ? sample + "_" + new_intervals.baseName : sample
            new_meta = [ id: id, sample: meta.sample, gender: meta.gender, status: meta.status, patient: meta.patient ]
            [new_meta, cram, crai, [], [], new_intervals]
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
            freebayes_vcf_gz = FREEBAYES.out.vcf
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
            freebayes_vcf_gz = CONCAT_VCF_FREEBAYES.out.vcf

            ch_versions = ch_versions.mix(BGZIP_FREEBAYES.out.versions)
            ch_versions = ch_versions.mix(CONCAT_VCF_FREEBAYES.out.versions)
        }

        freebayes_vcf = freebayes_vcf.mix(freebayes_vcf_gz)
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
            haplotypecaller_gvcf_gz = HAPLOTYPECALLER.out.vcf
            haplotypecaller_gvcf_gz_tbi = HAPLOTYPECALLER.out.tbi
            ch_versions = ch_versions.mix(TABIX_HAPLOTYPECALLER.out.versions)
        }else{
            BGZIP_HAPLOTYPECALLER(HAPLOTYPECALLER.out.vcf)

            BGZIP_HAPLOTYPECALLER.out.vcf.map{ meta, vcf ->
                new_meta = meta.clone()
                new_meta.id = new_meta.sample
                [new_meta, vcf]
            }.groupTuple(size: num_intervals)
            .set{haplotypecaller_gvcf_to_concat}

            CONCAT_VCF_HAPLOTYPECALLER(haplotypecaller_gvcf_to_concat, fasta_fai, intervals_bed_combine_gz)
            haplotypecaller_gvcf_gz = CONCAT_VCF_HAPLOTYPECALLER.out.vcf
            haplotypecaller_gvcf_gz_tbi = CONCAT_VCF_HAPLOTYPECALLER.out.tbi

            ch_versions = ch_versions.mix(BGZIP_HAPLOTYPECALLER.out.versions)
            ch_versions = ch_versions.mix(CONCAT_VCF_HAPLOTYPECALLER.out.versions)
        }

        haplotypecaller_gvcf_gz.join(haplotypecaller_gvcf_gz_tbi)
        .combine(intervals_bed_combine_gz_tbi)
        .map{
            meta, gvcf, gvf_tbi, intervals, intervals_tbi ->
            new_intervals = intervals.simpleName != "no_intervals" ? intervals : []
            new_intervals_tbi = intervals_tbi.simpleName != "no_intervals" ? intervals_tbi : []
            [meta, gvcf, gvf_tbi, new_intervals, new_intervals_tbi]
        }.set{haplotypecaller_gvcf_to_call}

        genotypegvcfs_vcf_gz = GENOTYPEGVCFS(
            haplotypecaller_gvcf_to_call,
            fasta,
            fasta_fai,
            dict,
            dbsnp,
            dbsnp_tbi)

        genotypegvcfs_vcf = genotypegvcfs_vcf.mix(genotypegvcfs_vcf_gz)

        if(joint_germline){
            run_haplotypecaller = false
            run_vqsr            = true //parameter?
            //some feedback from gavin
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

        haplotypecaller_gvcf = haplotypecaller_gvcf.mix(haplotypecaller_gvcf_gz)
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
            manta_candidate_small_indels_vcf = MANTA_GERMLINE.out.candidate_small_indels_vcf
            manta_candidate_sv_vcf           = MANTA_GERMLINE.out.candidate_sv_vcf
            manta_diploid_sv_vcf             = MANTA_GERMLINE.out.diploid_sv_vcf
        }else{

            BGZIP_MANTA_SV(MANTA_GERMLINE.out.candidate_small_indels_vcf)
            BGZIP_MANTA_SMALL_INDELS(MANTA_GERMLINE.out.candidate_sv_vcf)
            BGZIP_MANTA_DIPLOID(MANTA_GERMLINE.out.diploid_sv_vcf)

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

            BGZIP_MANTA_DIPLOID.out.vcf.map{ meta, vcf ->
                new_meta = meta.clone()
                new_meta.id = new_meta.sample
                [new_meta, vcf]
            }.groupTuple(size: num_intervals)
            .set{manta_diploid_vcf_to_concat}

            CONCAT_VCF_MANTA_SV(manta_sv_vcf_to_concat, fasta_fai, intervals_bed_combine_gz)
            CONCAT_VCF_MANTA_SMALL_INDELS(manta_small_indels_vcf_to_concat,fasta_fai, intervals_bed_combine_gz)
            CONCAT_VCF_MANTA_DIPLOID(manta_diploid_vcf_to_concat, fasta_fai, intervals_bed_combine_gz)

            manta_candidate_small_indels_vcf = CONCAT_VCF_MANTA_SV.out.vcf
            manta_candidate_sv_vcf           = CONCAT_VCF_MANTA_SMALL_INDELS.out.vcf
            manta_diploid_sv_vcf             = CONCAT_VCF_MANTA_DIPLOID.out.vcf

            ch_versions = ch_versions.mix(BGZIP_MANTA_SV.out.versions)
            ch_versions = ch_versions.mix(BGZIP_MANTA_SMALL_INDELS.out.versions)
            ch_versions = ch_versions.mix(BGZIP_MANTA_DIPLOID.out.versions)

            ch_versions = ch_versions.mix(CONCAT_VCF_MANTA_SV.out.versions)
            ch_versions = ch_versions.mix(CONCAT_VCF_MANTA_SMALL_INDELS.out.versions)
            ch_versions = ch_versions.mix(CONCAT_VCF_MANTA_DIPLOID.out.versions)

        }

        manta_vcf = manta_vcf.mix(manta_candidate_small_indels_vcf, manta_candidate_sv_vcf, manta_diploid_sv_vcf)
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
            strelka_vcf_gz = STRELKA_GERMLINE.out.vcf
            strelka_genome_vcf_gz = STRELKA_GERMLINE.out.genome_vcf

        }else{
            BGZIP_STRELKA(STRELKA_GERMLINE.out.vcf)
            BGZIP_STRELKA_GENOME(STRELKA_GERMLINE.out.genome_vcf)

            BGZIP_STRELKA.out.vcf.map{ meta, vcf ->
                new_meta = meta.clone()
                new_meta.id = new_meta.sample
                [new_meta, vcf]
            }.groupTuple(size: num_intervals)
            .set{strelka_vcf_to_concat}

            BGZIP_STRELKA_GENOME.out.vcf.map{ meta, vcf ->
                new_meta = meta.clone()
                new_meta.id = new_meta.sample
                [new_meta, vcf]
            }.groupTuple(size: num_intervals)
            .set{strelka_genome_vcf_to_concat}

            CONCAT_VCF_STRELKA(strelka_vcf_to_concat,fasta_fai, intervals_bed_combine_gz)
            CONCAT_VCF_STRELKA_GENOME(strelka_genome_vcf_to_concat,fasta_fai, intervals_bed_combine_gz)

            strelka_vcf_gz = CONCAT_VCF_STRELKA.out.vcf
            strelka_genome_vcf_gz = CONCAT_VCF_STRELKA_GENOME.out.vcf

            ch_versions = ch_versions.mix(BGZIP_STRELKA.out.versions)
            ch_versions = ch_versions.mix(CONCAT_VCF_STRELKA.out.versions)
        }

        strelka_vcf = strelka_vcf.mix(strelka_vcf_gz,strelka_genome_vcf_gz )
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
    deepvariant_vcf
    freebayes_vcf
    haplotypecaller_gvcf
    genotypegvcfs_vcf
    manta_vcf
    strelka_vcf

    versions = ch_versions
}
