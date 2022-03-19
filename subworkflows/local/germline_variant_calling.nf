//
// GERMLINE VARIANT CALLING
//

include { BGZIP as BGZIP_VC_FREEBAYES               } from '../../modules/local/bgzip'
include { BGZIP as BGZIP_VC_HAPLOTYPECALLER         } from '../../modules/local/bgzip'
include { BGZIP as BGZIP_VC_MANTA_DIPLOID           } from '../../modules/local/bgzip'
include { BGZIP as BGZIP_VC_MANTA_SMALL_INDELS      } from '../../modules/local/bgzip'
include { BGZIP as BGZIP_VC_MANTA_SV                } from '../../modules/local/bgzip'
include { BGZIP as BGZIP_VC_STRELKA                 } from '../../modules/local/bgzip'
include { BGZIP as BGZIP_VC_STRELKA_GENOME          } from '../../modules/local/bgzip'
include { CONCAT_VCF as CONCAT_FREEBAYES            } from '../../modules/local/concat_vcf/main'
include { CONCAT_VCF as CONCAT_HAPLOTYPECALLER      } from '../../modules/local/concat_vcf/main'
include { CONCAT_VCF as CONCAT_MANTA_DIPLOID        } from '../../modules/local/concat_vcf/main'
include { CONCAT_VCF as CONCAT_MANTA_SMALL_INDELS   } from '../../modules/local/concat_vcf/main'
include { CONCAT_VCF as CONCAT_MANTA_SV             } from '../../modules/local/concat_vcf/main'
include { CONCAT_VCF as CONCAT_STRELKA              } from '../../modules/local/concat_vcf/main'
include { CONCAT_VCF as CONCAT_STRELKA_GENOME       } from '../../modules/local/concat_vcf/main'
include { FREEBAYES                                 } from '../../modules/nf-core/modules/freebayes/main'
include { GATK4_GENOTYPEGVCFS as GENOTYPEGVCFS      } from '../../modules/nf-core/modules/gatk4/genotypegvcfs/main'
include { GATK4_HAPLOTYPECALLER as HAPLOTYPECALLER  } from '../../modules/nf-core/modules/gatk4/haplotypecaller/main'
include { GATK_JOINT_GERMLINE_VARIANT_CALLING       } from '../../subworkflows/nf-core/gatk4/joint_germline_variant_calling/main'
include { MANTA_GERMLINE                            } from '../../modules/local/manta/germline/main'
include { STRELKA_GERMLINE                          } from '../../modules/nf-core/modules/strelka/germline/main'
include { TABIX_BGZIPTABIX as TABIX_BGZIP_TIDDIT_SV } from '../../modules/nf-core/modules/tabix/bgziptabix/main'
include { TABIX_TABIX as TABIX_VC_FREEBAYES         } from '../../modules/nf-core/modules/tabix/tabix/main'
include { TABIX_TABIX as TABIX_VC_HAPLOTYPECALLER   } from '../../modules/nf-core/modules/tabix/tabix/main'
include { TIDDIT_SV                                 } from '../../modules/nf-core/modules/tiddit/sv/main'

include { RUN_DEEPVARIANT                           } from './variantcalling/deepvariant.nf'
workflow GERMLINE_VARIANT_CALLING {
    take:
        cram_recalibrated            // channel: [mandatory] cram
        dbsnp                        // channel: [mandatory] dbsnp
        dbsnp_tbi                    // channel: [mandatory] dbsnp_tbi
        dict                         // channel: [mandatory] dict
        fasta                        // channel: [mandatory] fasta
        fasta_fai                    // channel: [mandatory] fasta_fai
        intervals                    // channel: [mandatory] intervals/target regions
        intervals_bed_gz_tbi         // channel: [mandatory] intervals/target regions index zipped and indexed
        intervals_bed_combine_gz_tbi // channel: [mandatory] intervals/target regions index zipped and indexed in one file
        intervals_bed_combine_gz     // channel: [mandatory] intervals/target regions index zipped and indexed in one file
        num_intervals                // val: number of intervals that are used to parallelize exection, either based on capture kit or GATK recommended for WGS
        // joint_germline               // val: true/false on whether to run joint_germline calling, only works in combination with haplotypecaller at the moment

    main:

    ch_versions = Channel.empty()

    // Remap channel with intervals
    cram_recalibrated_intervals = cram_recalibrated.combine(intervals)
        .map{ meta, cram, crai, intervals ->
            sample = meta.sample
            new_intervals = intervals.baseName != "no_intervals" ? intervals : []
            id = new_intervals ? sample + "_" + new_intervals.baseName : sample
            [[ id: id, sample: meta.sample, gender: meta.gender, status: meta.status, patient: meta.patient ], cram, crai, new_intervals]
        }

    // Remap channel with gziped intervals + indexes
    cram_recalibrated_intervals_gz_tbi = cram_recalibrated.combine(intervals_bed_gz_tbi)
        .map{ meta, cram, crai, bed, tbi ->
            sample = meta.sample
            new_bed = bed.simpleName != "no_intervals" ? bed : []
            new_tbi = tbi.simpleName != "no_intervals" ? tbi : []
            id = new_bed ? sample + "_" + new_bed.simpleName : sample
            new_meta = [ id: id, sample: meta.sample, gender: meta.gender, status: meta.status, patient: meta.patient ]
            [new_meta, cram, crai, new_bed, new_tbi]
        }

    // DEEPVARIANT
    RUN_DEEPVARIANT(cram_recalibrated_intervals, fasta, fasta_fai, intervals_bed_combine_gz, num_intervals)
    // FREEBAYES

    // Remap channel for Freebayes
    cram_recalibrated_intervals_freebayes = cram_recalibrated.combine(intervals)
        .map{ meta, cram, crai, intervals ->
            sample = meta.sample
            new_intervals = intervals.baseName != "no_intervals" ? intervals : []
            id = new_intervals ? sample + "_" + new_intervals.baseName : sample
            new_meta = [ id: id, sample: meta.sample, gender: meta.gender, status: meta.status, patient: meta.patient ]
            [new_meta, cram, crai, [], [], new_intervals]
        }

    FREEBAYES(
        cram_recalibrated_intervals_freebayes,
        fasta,
        fasta_fai,
        [], [], [])

    // Only when no intervals
    TABIX_VC_FREEBAYES(FREEBAYES.out.vcf)

    // Only when using intervals
    BGZIP_VC_FREEBAYES(FREEBAYES.out.vcf)

    CONCAT_FREEBAYES(
        BGZIP_VC_FREEBAYES.out.vcf
            .map{ meta, vcf ->
                new_meta = meta.clone()
                new_meta.id = new_meta.sample
                [new_meta, vcf]
            }.groupTuple(size: num_intervals),
        fasta_fai,
        intervals_bed_combine_gz)

    freebayes_vcf = Channel.empty().mix(
        CONCAT_FREEBAYES.out.vcf,
        FREEBAYES.out.vcf.join(TABIX_VC_FREEBAYES.out.tbi))

    // HAPLOTYPECALLER

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

    // MANTA
    // TODO: Research if splitting by intervals is ok, we pretend for now it is fine.
    // Seems to be the consensus on upstream modules implementation too

    MANTA_GERMLINE(
        cram_recalibrated_intervals_gz_tbi,
        fasta,
        fasta_fai)

    // Figure out if using intervals or no_intervals
    MANTA_GERMLINE.out.candidate_small_indels_vcf.groupTuple(size: num_intervals)
        .branch{
            intervals:    it[1].size() > 1
            no_intervals: it[1].size() == 1
        }.set{manta_small_indels_vcf}

    MANTA_GERMLINE.out.candidate_sv_vcf.groupTuple(size: num_intervals)
        .branch{
            intervals:    it[1].size() > 1
            no_intervals: it[1].size() == 1
        }.set{manta_sv_vcf}

    MANTA_GERMLINE.out.diploid_sv_vcf.groupTuple(size: num_intervals)
        .branch{
            intervals:    it[1].size() > 1
            no_intervals: it[1].size() == 1
        }.set{manta_diploid_sv_vcf}

    // Only when using intervals
    BGZIP_VC_MANTA_DIPLOID(MANTA_GERMLINE.out.diploid_sv_vcf)

    CONCAT_MANTA_DIPLOID(
        BGZIP_VC_MANTA_DIPLOID.out.vcf
            .map{ meta, vcf ->
                new_meta = meta.clone()
                new_meta.id = new_meta.sample
                [new_meta, vcf]
            }.groupTuple(size: num_intervals),
        fasta_fai,
        intervals_bed_combine_gz)

    BGZIP_VC_MANTA_SMALL_INDELS(MANTA_GERMLINE.out.candidate_small_indels_vcf)

    CONCAT_MANTA_SMALL_INDELS(
        BGZIP_VC_MANTA_SMALL_INDELS.out.vcf
            .map{ meta, vcf ->
                new_meta = meta.clone()
                new_meta.id = new_meta.sample
                [new_meta, vcf]
            }.groupTuple(size: num_intervals),
        fasta_fai,
        intervals_bed_combine_gz)

    BGZIP_VC_MANTA_SV(MANTA_GERMLINE.out.candidate_sv_vcf)

    CONCAT_MANTA_SV(
        BGZIP_VC_MANTA_SV.out.vcf
            .map{ meta, vcf ->
                new_meta = meta.clone()
                new_meta.id = new_meta.sample
                [new_meta, vcf]
            }.groupTuple(size: num_intervals),
        fasta_fai,
        intervals_bed_combine_gz)

    manta_vcf = Channel.empty().mix(
        CONCAT_MANTA_DIPLOID.out.vcf,
        CONCAT_MANTA_SMALL_INDELS.out.vcf,
        CONCAT_MANTA_SV.out.vcf,
        manta_diploid_sv_vcf.no_intervals,
        manta_small_indels_vcf.no_intervals,
        manta_sv_vcf.no_intervals)

    // STRELKA
    // TODO: Research if splitting by intervals is ok, we pretend for now it is fine.
    // Seems to be the consensus on upstream modules implementation too

    STRELKA_GERMLINE(
        cram_recalibrated_intervals_gz_tbi,
        fasta,
        fasta_fai)

    // Figure out if using intervals or no_intervals
    STRELKA_GERMLINE.out.vcf.groupTuple(size: num_intervals)
        .branch{
            intervals:    it[1].size() > 1
            no_intervals: it[1].size() == 1
        }.set{strelka_vcf}

    STRELKA_GERMLINE.out.genome_vcf.groupTuple(size: num_intervals)
        .branch{
            intervals:    it[1].size() > 1
            no_intervals: it[1].size() == 1
        }.set{strelka_genome_vcf}

    // Only when using intervals
    BGZIP_VC_STRELKA(STRELKA_GERMLINE.out.vcf)

    CONCAT_STRELKA(
        BGZIP_VC_STRELKA.out.vcf
            .map{ meta, vcf ->
                new_meta = meta.clone()
                new_meta.id = new_meta.sample
                [new_meta, vcf]
            }.groupTuple(size: num_intervals),
        fasta_fai,
        intervals_bed_combine_gz)

    BGZIP_VC_STRELKA_GENOME(STRELKA_GERMLINE.out.genome_vcf)

    CONCAT_STRELKA_GENOME(
        BGZIP_VC_STRELKA_GENOME.out.vcf
            .map{ meta, vcf ->
                new_meta = meta.clone()
                new_meta.id = new_meta.sample
                [new_meta, vcf]
            }.groupTuple(size: num_intervals),
        fasta_fai,
        intervals_bed_combine_gz)

    strelka_vcf = Channel.empty().mix(
        CONCAT_STRELKA.out.vcf,
        CONCAT_STRELKA_GENOME.out.vcf,
        strelka_genome_vcf.no_intervals,
        strelka_vcf.no_intervals)

    // if (tools.contains('tiddit')) {
    //     TODO: Update tiddit on bioconda, the current version does not support cram usage, needs newest version:
    //     https://github.com/SciLifeLab/TIDDIT/issues/82#issuecomment-1022103264
    //     Issue opened, either this week or end of february

    //     TIDDIT_SV(
    //         cram_recalibrated,
    //         fasta,
    //         fasta_fai
    //     )

    //     TABIX_BGZIP_TIDDIT_SV(TIDDIT_SV.out.vcf)
    //     tiddit_vcf_gz_tbi = TABIX_BGZIP_TIDDIT_SV.out.gz_tbi
    //     tiddit_ploidy     = TIDDIT_SV.out.ploidy
    //     tiddit_signals    = TIDDIT_SV.out.signals
    //     tiddit_wig        = TIDDIT_SV.out.wig
    //     tiddit_gc_wig     = TIDDIT_SV.out.gc_wig

    //     ch_versions = ch_versions.mix(TABIX_BGZIP_TIDDIT_SV.out.versions)
    //     ch_versions = ch_versions.mix(TIDDIT_SV.out.versions)
    // }


    ch_versions = ch_versions.mix(BGZIP_VC_FREEBAYES.out.versions)
    ch_versions = ch_versions.mix(BGZIP_VC_HAPLOTYPECALLER.out.versions)
    ch_versions = ch_versions.mix(BGZIP_VC_MANTA_DIPLOID.out.versions)
    ch_versions = ch_versions.mix(BGZIP_VC_MANTA_SMALL_INDELS.out.versions)
    ch_versions = ch_versions.mix(BGZIP_VC_MANTA_SV.out.versions)
    ch_versions = ch_versions.mix(BGZIP_VC_STRELKA.out.versions)
    ch_versions = ch_versions.mix(CONCAT_FREEBAYES.out.versions)
    ch_versions = ch_versions.mix(CONCAT_HAPLOTYPECALLER.out.versions)
    ch_versions = ch_versions.mix(CONCAT_MANTA_DIPLOID.out.versions)
    ch_versions = ch_versions.mix(CONCAT_MANTA_SMALL_INDELS.out.versions)
    ch_versions = ch_versions.mix(CONCAT_MANTA_SV.out.versions)
    ch_versions = ch_versions.mix(CONCAT_STRELKA.out.versions)
    ch_versions = ch_versions.mix(FREEBAYES.out.versions)
    ch_versions = ch_versions.mix(GENOTYPEGVCFS.out.versions)
    ch_versions = ch_versions.mix(HAPLOTYPECALLER.out.versions)
    ch_versions = ch_versions.mix(MANTA_GERMLINE.out.versions)
    ch_versions = ch_versions.mix(STRELKA_GERMLINE.out.versions)
    ch_versions = ch_versions.mix(TABIX_VC_FREEBAYES.out.versions)
    ch_versions = ch_versions.mix(TABIX_VC_HAPLOTYPECALLER.out.versions)

    ch_versions = ch_versions.mix(RUN_DEEPVARIANT.out.versions)

    emit:
    deepvariant_vcf = RUN_DEEPVARIANT.out.deepvariant_vcf
    freebayes_vcf
    haplotypecaller_gvcf
    genotype_gvcf
    manta_vcf
    strelka_vcf

    versions = ch_versions
}
