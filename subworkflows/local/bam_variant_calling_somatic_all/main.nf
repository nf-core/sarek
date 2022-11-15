//
// PAIRED VARIANT CALLING
//

include { BAM_VARIANT_CALLING_CNVKIT                    } from '../bam_variant_calling_cnvkit/main'
include { BAM_VARIANT_CALLING_FREEBAYES                 } from '../bam_variant_calling_freebayes/main'
include { BAM_VARIANT_CALLING_MPILEUP as MPILEUP_NORMAL } from '../bam_variant_calling_mpileup/main'
include { BAM_VARIANT_CALLING_MPILEUP as MPILEUP_TUMOR  } from '../bam_variant_calling_mpileup/main'
include { BAM_VARIANT_CALLING_SOMATIC_ASCAT             } from '../bam_variant_calling_somatic_ascat/main'
include { BAM_VARIANT_CALLING_SOMATIC_CONTROLFREEC      } from '../bam_variant_calling_somatic_controlfreec/main'
include { BAM_VARIANT_CALLING_SOMATIC_MANTA             } from '../bam_variant_calling_somatic_manta/main'
include { BAM_VARIANT_CALLING_SOMATIC_MUTECT2           } from '../bam_variant_calling_somatic_mutect2/main'
include { BAM_VARIANT_CALLING_SOMATIC_STRELKA           } from '../bam_variant_calling_somatic_strelka/main'
include { BAM_VARIANT_CALLING_SOMATIC_TIDDIT            } from '../bam_variant_calling_somatic_tiddit/main'
include { MSISENSORPRO_MSI_SOMATIC                      } from '../../../modules/nf-core/msisensorpro/msi_somatic/main'

workflow BAM_VARIANT_CALLING_SOMATIC_ALL {
    take:
        tools                         // Mandatory, list of tools to apply
        cram_pair                     // channel: [mandatory] cram
        bwa                           // channel: [optional] bwa
        cf_chrom_len                  // channel: [optional] controlfreec length file
        chr_files
        dbsnp                         // channel: [mandatory] dbsnp
        dbsnp_tbi                     // channel: [mandatory] dbsnp_tbi
        dict                          // channel: [mandatory] dict
        fasta                         // channel: [mandatory] fasta
        fasta_fai                     // channel: [mandatory] fasta_fai
        germline_resource             // channel: [optional]  germline_resource
        germline_resource_tbi         // channel: [optional]  germline_resource_tbi
        intervals                     // channel: [mandatory] intervals/target regions
        intervals_bed_gz_tbi          // channel: [mandatory] intervals/target regions index zipped and indexed
        intervals_bed_combined        // channel: [mandatory] intervals/target regions in one file unzipped
        mappability
        msisensorpro_scan             // channel: [optional]  msisensorpro_scan
        panel_of_normals              // channel: [optional]  panel_of_normals
        panel_of_normals_tbi          // channel: [optional]  panel_of_normals_tbi
        allele_files                  // channel: [optional]  ascat allele files
        loci_files                    // channel: [optional]  ascat loci files
        gc_file                       // channel: [optional]  ascat gc content file
        rt_file                       // channel: [optional]  ascat rt file

    main:

    ch_versions          = Channel.empty()

    //TODO: Temporary until the if's can be removed and printing to terminal is prevented with "when" in the modules.config
    freebayes_vcf        = Channel.empty()
    manta_vcf            = Channel.empty()
    strelka_vcf          = Channel.empty()
    msisensorpro_output  = Channel.empty()
    mutect2_vcf          = Channel.empty()
    tiddit_vcf           = Channel.empty()

    // Remap channel with intervals
    cram_pair_intervals = cram_pair.combine(intervals)
        .map{ meta, normal_cram, normal_crai, tumor_cram, tumor_crai, intervals, num_intervals ->
            //If no interval file provided (0) then add empty list
            intervals_new = num_intervals == 0 ? [] : intervals

            [[
                id:             meta.tumor_id + "_vs_" + meta.normal_id,
                normal_id:      meta.normal_id,
                num_intervals:  num_intervals,
                patient:        meta.patient,
                sex:            meta.sex,
                tumor_id:       meta.tumor_id,
            ],
            normal_cram, normal_crai, tumor_cram, tumor_crai, intervals_new]
        }

    // Remap channel with gzipped intervals + indexes
    cram_pair_intervals_gz_tbi = cram_pair.combine(intervals_bed_gz_tbi)
        .map{ meta, normal_cram, normal_crai, tumor_cram, tumor_crai, bed_tbi, num_intervals ->

            //If no interval file provided (0) then add empty list
            bed_new = num_intervals == 0 ? [] : bed_tbi[0]
            tbi_new = num_intervals == 0 ? [] : bed_tbi[1]

            [[
                id:             meta.tumor_id + "_vs_" + meta.normal_id,
                normal_id:      meta.normal_id,
                num_intervals:  num_intervals,
                patient:        meta.patient,
                sex:            meta.sex,
                tumor_id:       meta.tumor_id,
            ],
            normal_cram, normal_crai, tumor_cram, tumor_crai, bed_new, tbi_new]

        }

    if (tools.split(',').contains('ascat')){

        BAM_VARIANT_CALLING_SOMATIC_ASCAT(
            cram_pair,
            allele_files,
            loci_files,
            intervals_bed_combined,
            fasta,
            gc_file,
            rt_file
        )

        ch_versions = ch_versions.mix(BAM_VARIANT_CALLING_SOMATIC_ASCAT.out.versions)

    }

    if (tools.split(',').contains('controlfreec')){
        cram_normal_intervals_no_index = cram_pair_intervals
                    .map {meta, normal_cram, normal_crai, tumor_cram, tumor_crai, intervals ->
                            [meta, normal_cram, intervals]
                        }

        cram_tumor_intervals_no_index = cram_pair_intervals
                    .map {meta, normal_cram, normal_crai, tumor_cram, tumor_crai, intervals ->
                            [meta, tumor_cram, intervals]
                        }

        MPILEUP_NORMAL(
            cram_normal_intervals_no_index,
            fasta
        )

        MPILEUP_TUMOR(
            cram_tumor_intervals_no_index,
            fasta
        )

        mpileup_normal = MPILEUP_NORMAL.out.mpileup
        mpileup_tumor = MPILEUP_TUMOR.out.mpileup

        controlfreec_input = mpileup_normal.cross(mpileup_tumor)
        .map{ normal, tumor ->
            [normal[0], normal[1], tumor[1], [], [], [], []]
        }

        length_file = cf_chrom_len ?: fasta_fai
        BAM_VARIANT_CALLING_SOMATIC_CONTROLFREEC(
            controlfreec_input,
            fasta,
            length_file,
            dbsnp,
            dbsnp_tbi,
            chr_files,
            mappability,
            intervals_bed_combined
        )

        ch_versions = ch_versions.mix(MPILEUP_NORMAL.out.versions)
        ch_versions = ch_versions.mix(MPILEUP_TUMOR.out.versions)
        ch_versions = ch_versions.mix(BAM_VARIANT_CALLING_SOMATIC_CONTROLFREEC.out.versions)
    }

    if (tools.split(',').contains('cnvkit')){
        cram_pair_cnvkit_somatic = cram_pair
            .map{meta, normal_cram, normal_crai, tumor_cram, tumor_crai ->
                [meta, tumor_cram, normal_cram]
            }

        BAM_VARIANT_CALLING_CNVKIT(
            cram_pair_cnvkit_somatic,
            fasta,
            fasta_fai,
            intervals_bed_combined,
            []
        )

        ch_versions = ch_versions.mix(BAM_VARIANT_CALLING_CNVKIT.out.versions)
    }

    if (tools.split(',').contains('freebayes')){

        BAM_VARIANT_CALLING_FREEBAYES(
            cram_pair_intervals,
            dict,
            fasta,
            fasta_fai
        )

        freebayes_vcf = BAM_VARIANT_CALLING_FREEBAYES.out.freebayes_vcf
        ch_versions   = ch_versions.mix(BAM_VARIANT_CALLING_FREEBAYES.out.versions)
    }

    if (tools.split(',').contains('manta')) {
        BAM_VARIANT_CALLING_SOMATIC_MANTA(
            cram_pair_intervals_gz_tbi,
            dict,
            fasta,
            fasta_fai
        )

        manta_vcf                            = BAM_VARIANT_CALLING_SOMATIC_MANTA.out.manta_vcf
        manta_candidate_small_indels_vcf     = BAM_VARIANT_CALLING_SOMATIC_MANTA.out.manta_candidate_small_indels_vcf
        manta_candidate_small_indels_vcf_tbi = BAM_VARIANT_CALLING_SOMATIC_MANTA.out.manta_candidate_small_indels_vcf_tbi
        ch_versions                          = ch_versions.mix(BAM_VARIANT_CALLING_SOMATIC_MANTA.out.versions)
    }

    if (tools.split(',').contains('strelka')) {

        if (tools.split(',').contains('manta')) {
            cram_pair_strelka = cram_pair.join(manta_candidate_small_indels_vcf)
                                        .join(manta_candidate_small_indels_vcf_tbi)
                                        .combine(intervals_bed_gz_tbi)
                                        .map{
                                            meta, normal_cram, normal_crai, tumor_cram, tumor_crai, vcf, vcf_tbi, bed_tbi, num_intervals ->
                                            //If no interval file provided (0) then add empty list
                                            bed_new = num_intervals == 0 ? [] : bed_tbi[0]
                                            tbi_new = num_intervals == 0 ? [] : bed_tbi[1]

                                            [[
                                                id:             meta.tumor_id + "_vs_" + meta.normal_id,
                                                normal_id:      meta.normal_id,
                                                num_intervals:  num_intervals,
                                                patient:        meta.patient,
                                                sex:            meta.sex,
                                                tumor_id:       meta.tumor_id,
                                            ],
                                            normal_cram, normal_crai, tumor_cram, tumor_crai, vcf, vcf_tbi, bed_new, tbi_new]
                                        }

        } else {
            cram_pair_strelka = cram_pair_intervals_gz_tbi.map{
                    meta, normal_cram, normal_crai, tumor_cram, tumor_crai, bed, tbi ->
                    [meta, normal_cram, normal_crai, tumor_cram, tumor_crai, [], [], bed, tbi]
            }
        }

        BAM_VARIANT_CALLING_SOMATIC_STRELKA(
            cram_pair_strelka,
            dict,
            fasta,
            fasta_fai
        )

        strelka_vcf = Channel.empty().mix(BAM_VARIANT_CALLING_SOMATIC_STRELKA.out.strelka_vcf)
        ch_versions = ch_versions.mix(BAM_VARIANT_CALLING_SOMATIC_STRELKA.out.versions)
    }

    if (tools.split(',').contains('msisensorpro')) {

        cram_pair_msisensor = cram_pair.combine(intervals_bed_combined)
        MSISENSORPRO_MSI_SOMATIC(cram_pair_msisensor, fasta, msisensorpro_scan)
        ch_versions = ch_versions.mix(MSISENSORPRO_MSI_SOMATIC.out.versions)
        msisensorpro_output = msisensorpro_output.mix(MSISENSORPRO_MSI_SOMATIC.out.output_report)
    }

    if (tools.split(',').contains('mutect2')) {
        cram_pair_mutect2 = cram_pair_intervals.map{ meta, normal_cram, normal_crai, tumor_cram, tumor_crai, intervals ->
                                [meta, [normal_cram, tumor_cram], [normal_crai, tumor_crai], intervals]
                            }

        BAM_VARIANT_CALLING_SOMATIC_MUTECT2(
            cram_pair_mutect2,
            fasta,
            fasta_fai,
            dict,
            germline_resource,
            germline_resource_tbi,
            panel_of_normals,
            panel_of_normals_tbi
        )

        mutect2_vcf = BAM_VARIANT_CALLING_SOMATIC_MUTECT2.out.filtered_vcf
        ch_versions = ch_versions.mix(BAM_VARIANT_CALLING_SOMATIC_MUTECT2.out.versions)
    }

    //TIDDIT
    if (tools.split(',').contains('tiddit')){
        cram_normal = cram_pair.map{meta, normal_cram, normal_crai, tumor_cram, tumor_crai ->
            [meta, normal_cram, normal_crai]
        }
        cram_tumor = cram_pair.map{meta, normal_cram, normal_crai, tumor_cram, tumor_crai ->
            [meta, tumor_cram, tumor_crai]
        }

        BAM_VARIANT_CALLING_SOMATIC_TIDDIT(cram_normal, cram_tumor, fasta.map{ it -> [[id:it[0].baseName], it] }, bwa)
        tiddit_vcf = BAM_VARIANT_CALLING_SOMATIC_TIDDIT.out.tiddit_vcf
        ch_versions = ch_versions.mix(BAM_VARIANT_CALLING_SOMATIC_TIDDIT.out.versions)
    }

    emit:
    freebayes_vcf
    manta_vcf
    msisensorpro_output
    mutect2_vcf
    strelka_vcf
    tiddit_vcf

    versions    = ch_versions
}
