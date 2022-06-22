//
// PAIRED VARIANT CALLING
//
include { GATK_TUMOR_NORMAL_SOMATIC_VARIANT_CALLING } from '../../subworkflows/nf-core/gatk4/tumor_normal_somatic_variant_calling/main'
include { MSISENSORPRO_MSI_SOMATIC                  } from '../../modules/nf-core/modules/msisensorpro/msi_somatic/main'
include { RUN_CONTROLFREEC_SOMATIC                  } from '../nf-core/variantcalling/controlfreec/somatic/main.nf'
include { RUN_FREEBAYES as RUN_FREEBAYES_SOMATIC    } from '../nf-core/variantcalling/freebayes/main.nf'
include { RUN_MANTA_SOMATIC                         } from '../nf-core/variantcalling/manta/somatic/main.nf'
include { RUN_STRELKA_SOMATIC                       } from '../nf-core/variantcalling/strelka/somatic/main.nf'
include { RUN_CNVKIT_SOMATIC                        } from '../nf-core/variantcalling/cnvkit/somatic/main.nf'
include { RUN_MPILEUP as RUN_MPILEUP_NORMAL         } from '../nf-core/variantcalling/mpileup/main'
include { RUN_MPILEUP as RUN_MPILEUP_TUMOR          } from '../nf-core/variantcalling/mpileup/main'
include { RUN_ASCAT_SOMATIC                         } from '../nf-core/variantcalling/ascat/main'

workflow PAIR_VARIANT_CALLING {
    take:
        tools                         // Mandatory, list of tools to apply
        cram_pair                     // channel: [mandatory] cram
        dbsnp                         // channel: [mandatory] dbsnp
        dbsnp_tbi                     // channel: [mandatory] dbsnp_tbi
        dict                          // channel: [mandatory] dict
        fasta                         // channel: [mandatory] fasta
        fasta_fai                     // channel: [mandatory] fasta_fai
        intervals                     // channel: [mandatory] intervals/target regions
        intervals_bed_gz_tbi          // channel: [mandatory] intervals/target regions index zipped and indexed
        intervals_bed_combined        // channel: [mandatory] intervals/target regions in one file unzipped
        msisensorpro_scan             // channel: [optional]  msisensorpro_scan
        germline_resource             // channel: [optional]  germline_resource
        germline_resource_tbi         // channel: [optional]  germline_resource_tbi
        panel_of_normals              // channel: [optional]  panel_of_normals
        panel_of_normals_tbi          // channel: [optional]  panel_of_normals_tbi
        chr_files
        mappability
        allele_files                  // channel: [optional]  ascat allele files
        loci_files                    // channel: [optional]  ascat loci files
        gc_file                       // channel: [optional]  ascat gc content file
        rt_file                       // channel: [optional]  ascat rt file

    main:

    ch_versions          = Channel.empty()
    cram_pair.view()

    //TODO: Temporary until the if's can be removed and printing to terminal is prevented with "when" in the modules.config
    freebayes_vcf        = Channel.empty()
    manta_vcf            = Channel.empty()
    strelka_vcf          = Channel.empty()
    msisensorpro_output  = Channel.empty()
    mutect2_vcf          = Channel.empty()

    // Remap channel with intervals
    cram_pair_intervals = cram_pair.combine(intervals)
        .map{ meta, normal_cram, normal_crai, tumor_cram, tumor_crai, intervals, num_intervals ->
            //If no interval file provided (0) then add empty list
            intervals_new = num_intervals == 0 ? [] : intervals

            [[patient:meta.patient, normal_id:meta.normal_id, tumor_id:meta.tumor_id, gender:meta.gender, id: meta.tumor_id + "_vs_" + meta.normal_id, num_intervals:num_intervals],
            normal_cram, normal_crai, tumor_cram, tumor_crai, intervals_new]
        }

    // Remap channel with gzipped intervals + indexes
    cram_pair_intervals_gz_tbi = cram_pair.combine(intervals_bed_gz_tbi)
        .map{ meta, normal_cram, normal_crai, tumor_cram, tumor_crai, bed_tbi, num_intervals ->

            //If no interval file provided (0) then add empty list
            bed_new = num_intervals == 0 ? [] : bed_tbi[0]
            tbi_new = num_intervals == 0 ? [] : bed_tbi[1]

            [[patient:meta.patient, normal_id:meta.normal_id, tumor_id:meta.tumor_id, gender:meta.gender, id: meta.tumor_id + "_vs_" + meta.normal_id, num_intervals:num_intervals],
            normal_cram, normal_crai, tumor_cram, tumor_crai, bed_new, tbi_new]

        }

    if (tools.contains('ascat')){

        cram_pair.view()

        RUN_ASCAT_SOMATIC(cram_pair,
                        allele_files,
                        loci_files,
                        gc_file,
                        rt_file,
                        [],
                        fasta)

        ch_versions = ch_versions.mix(RUN_ASCAT_SOMATIC.out.versions)

    }

    if (tools.contains('controlfreec')){
        cram_normal_intervals_no_index = cram_pair_intervals
                    .map {meta, normal_cram, normal_crai, tumor_cram, tumor_crai, intervals ->
                            [meta, normal_cram, intervals]
                        }

        cram_tumor_intervals_no_index = cram_pair_intervals
                    .map {meta, normal_cram, normal_crai, tumor_cram, tumor_crai, intervals ->
                            [meta, tumor_cram, intervals]
                        }
        RUN_MPILEUP_NORMAL(cram_normal_intervals_no_index, fasta)
        mpileup_normal = RUN_MPILEUP_NORMAL.out.mpileup
        RUN_MPILEUP_TUMOR(cram_tumor_intervals_no_index, fasta)
        mpileup_tumor = RUN_MPILEUP_TUMOR.out.mpileup
        ch_versions = ch_versions.mix(RUN_MPILEUP_NORMAL.out.versions)
        ch_versions = ch_versions.mix(RUN_MPILEUP_TUMOR.out.versions)

        controlfreec_input = mpileup_normal.cross(mpileup_tumor)
        .map{ normal, tumor ->
            [normal[0], normal[1], tumor[1], [], [], [], []]
        }
        RUN_CONTROLFREEC_SOMATIC(controlfreec_input,
                        fasta,
                        fasta_fai,
                        dbsnp,
                        dbsnp_tbi,
                        chr_files,
                        mappability,
                        intervals_bed_combined)
        ch_versions = ch_versions.mix(RUN_CONTROLFREEC_SOMATIC.out.versions)
    }

    if (tools.contains('cnvkit')){
        cram_pair_cnvkit_somatic = cram_pair
            .map{meta, normal_cram, normal_crai, tumor_cram, tumor_crai ->
                [meta, tumor_cram, normal_cram]
            }

        RUN_CNVKIT_SOMATIC( cram_pair_cnvkit_somatic,
                            fasta,
                            fasta_fai,
                            intervals_bed_combined,
                            [])
    }

    if (tools.contains('freebayes')){
        RUN_FREEBAYES_SOMATIC(cram_pair_intervals, dict, fasta, fasta_fai)

        freebayes_vcf = RUN_FREEBAYES_SOMATIC.out.freebayes_vcf
        ch_versions   = ch_versions.mix(RUN_FREEBAYES_SOMATIC.out.versions)
    }

    if (tools.contains('manta')) {
        RUN_MANTA_SOMATIC(  cram_pair_intervals_gz_tbi,
                            dict,
                            fasta,
                            fasta_fai)

        manta_vcf                            = RUN_MANTA_SOMATIC.out.manta_vcf
        manta_candidate_small_indels_vcf     = RUN_MANTA_SOMATIC.out.manta_candidate_small_indels_vcf
        manta_candidate_small_indels_vcf_tbi = RUN_MANTA_SOMATIC.out.manta_candidate_small_indels_vcf_tbi
        ch_versions                          = ch_versions.mix(RUN_MANTA_SOMATIC.out.versions)
    }

    if (tools.contains('strelka')) {

        if (tools.contains('manta')) {
            cram_pair_strelka = cram_pair.join(manta_candidate_small_indels_vcf)
                                        .join(manta_candidate_small_indels_vcf_tbi)
                                        .combine(intervals_bed_gz_tbi)
                                        .map{
                                            meta, normal_cram, normal_crai, tumor_cram, tumor_crai, vcf, vcf_tbi, bed_tbi, num_intervals ->
                                            //If no interval file provided (0) then add empty list
                                            bed_new = num_intervals == 0 ? [] : bed_tbi[0]
                                            tbi_new = num_intervals == 0 ? [] : bed_tbi[1]

                                            [[patient:meta.patient, normal_id:meta.normal_id, tumor_id:meta.tumor_id, gender:meta.gender, id:meta.tumor_id + "_vs_" + meta.normal_id, num_intervals:num_intervals],
                                            normal_cram, normal_crai, tumor_cram, tumor_crai, vcf, vcf_tbi, bed_new, tbi_new]
                                        }

        } else {
            cram_pair_strelka = cram_pair_intervals_gz_tbi.map{
                    meta, normal_cram, normal_crai, tumor_cram, tumor_crai, bed, tbi ->
                    [meta, normal_cram, normal_crai, tumor_cram, tumor_crai, [], [], bed, tbi]
            }
        }

        RUN_STRELKA_SOMATIC(cram_pair_strelka,
                            dict,
                            fasta,
                            fasta_fai)

        strelka_vcf  = Channel.empty().mix(RUN_STRELKA_SOMATIC.out.strelka_vcf)
        ch_versions = ch_versions.mix(RUN_STRELKA_SOMATIC.out.versions)
    }

    if (tools.contains('msisensorpro')) {

        cram_pair_msisensor = cram_pair.combine(intervals_bed_combined)
        MSISENSORPRO_MSI_SOMATIC(cram_pair_msisensor, fasta, msisensorpro_scan)
        ch_versions = ch_versions.mix(MSISENSORPRO_MSI_SOMATIC.out.versions)
        msisensorpro_output = msisensorpro_output.mix(MSISENSORPRO_MSI_SOMATIC.out.output_report)
    }

    if (tools.contains('mutect2')) {
        cram_pair_mutect2 = cram_pair_intervals.map{ meta, normal_cram, normal_crai, tumor_cram, tumor_crai, intervals ->
                                [meta, [normal_cram, tumor_cram], [normal_crai, tumor_crai], intervals]
                            }

        GATK_TUMOR_NORMAL_SOMATIC_VARIANT_CALLING(
            cram_pair_mutect2,
            fasta,
            fasta_fai,
            dict,
            germline_resource,
            germline_resource_tbi,
            panel_of_normals,
            panel_of_normals_tbi)

        mutect2_vcf = GATK_TUMOR_NORMAL_SOMATIC_VARIANT_CALLING.out.filtered_vcf
        ch_versions = ch_versions.mix(GATK_TUMOR_NORMAL_SOMATIC_VARIANT_CALLING.out.versions)
    }

    // if (tools.contains('tiddit')) {
    // }

    emit:
    freebayes_vcf
    manta_vcf
    msisensorpro_output
    mutect2_vcf
    strelka_vcf
    versions    = ch_versions
}
