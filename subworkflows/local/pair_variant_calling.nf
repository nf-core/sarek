//
// PAIRED VARIANT CALLING
//
include { GATK_TUMOR_NORMAL_SOMATIC_VARIANT_CALLING } from '../../subworkflows/nf-core/gatk4/tumor_normal_somatic_variant_calling/main'
include { MSISENSORPRO_MSI_SOMATIC                  } from '../../modules/nf-core/modules/msisensorpro/msi_somatic/main'
include { RUN_MANTA_SOMATIC                         } from './variantcalling/manta_somatic.nf'
include { RUN_STRELKA_SOMATIC                       } from './variantcalling/strelka_somatic.nf'

workflow PAIR_VARIANT_CALLING {
    take:
        tools
        cram_pair                     // channel: [mandatory] cram
        dbsnp                         // channel: [mandatory] dbsnp
        dbsnp_tbi                     // channel: [mandatory] dbsnp_tbi
        dict                          // channel: [mandatory] dict
        fasta                         // channel: [mandatory] fasta
        fasta_fai                     // channel: [mandatory] fasta_fai
        intervals                     // channel: [mandatory] intervals/target regions
        intervals_bed_gz_tbi          // channel: [mandatory] intervals/target regions index zipped and indexed
        intervals_bed_combined_gz_tbi // channel: [mandatory] intervals/target regions index zipped and indexed
        intervals_bed_combine_gz      // channel: [mandatory] intervals/target regions index zipped and indexed in one file
        num_intervals                 // val: number of intervals that are used to parallelize exection, either based on capture kit or GATK recommended for WGS
        no_intervals
        msisensorpro_scan             // channel: [optional]  msisensorpro_scan
        germline_resource             // channel: [optional]  germline_resource
        germline_resource_tbi         // channel: [optional]  germline_resource_tbi
        panel_of_normals              // channel: [optional]  panel_of_normals
        panel_of_normals_tbi          // channel: [optional]  panel_of_normals_tbi

    main:

    ch_versions          = Channel.empty()

    //TODO: Temporary until the if's can be removed and printing to terminal is prevented with "when" in the modules.config
    manta_vcf            = Channel.empty()
    strelka_vcf          = Channel.empty()
    msisensorpro_output  = Channel.empty()
    mutect2_vcf          = Channel.empty()

    cram_pair_intervals_gz_tbi = cram_pair.combine(intervals_bed_gz_tbi)
        .map{ meta, normal_cram, normal_crai, tumor_cram, tumor_crai, bed, tbi ->
            normal_id = meta.normal_id
            tumor_id = meta.tumor_id

            new_bed = bed.simpleName != "no_intervals" ? bed : []
            new_tbi = tbi.simpleName != "no_intervals" ? tbi : []
            id = bed.simpleName != "no_intervals" ? tumor_id + "_vs_" + normal_id + "_" + bed.simpleName : tumor_id + "_vs_" + normal_id
            new_meta = [ id: id, normal_id: meta.normal_id, tumor_id: meta.tumor_id, gender: meta.gender, patient: meta.patient]
            [new_meta, normal_cram, normal_crai, tumor_cram, tumor_crai, new_bed, new_tbi]
        }

    cram_pair_intervals = cram_pair.combine(intervals)
        .map{ meta, normal_cram, normal_crai, tumor_cram, tumor_crai, intervals ->
            normal_id = meta.normal_id
            tumor_id = meta.tumor_id
            new_intervals = intervals.baseName != "no_intervals" ? intervals : []
            id = new_intervals ? tumor_id + "_vs_" + normal_id + "_" + new_intervals.baseName : tumor_id + "_vs_" + normal_id
            new_meta = [ id: id, normal_id: meta.normal_id, tumor_id: meta.tumor_id, gender: meta.gender, patient: meta.patient ]
            [new_meta, normal_cram, normal_crai, tumor_cram, tumor_crai, intervals]
        }

    if (tools.contains('manta')) {
        RUN_MANTA_SOMATIC(  cram_pair_intervals_gz_tbi,
                            fasta,
                            fasta_fai,
                            intervals_bed_combine_gz,
                            num_intervals)

        manta_vcf   = RUN_MANTA_SOMATIC.out.manta_vcf
        ch_versions = ch_versions.mix(RUN_MANTA_SOMATIC.out.versions)
    }

    if (tools.contains('strelka')) {

        if (tools.contains('manta')) {
            cram_pair_strelka = intervals_bed_gz_tbi.join(manta_somatic_sv_vcf).map{
                meta, normal_cram, normal_crai, tumor_cram, tumor_crai, bed, tbi, manta_vcf, manta_tbi ->
                [meta, normal_cram, normal_crai, tumor_cram, tumor_crai, manta_vcf, manta_tbi, bed, tbi]
            }
        } else {
            cram_pair_strelka = cram_pair_intervals_gz_tbi.map{
                    meta, normal_cram, normal_crai, tumor_cram, tumor_crai, bed, tbi ->
                    [meta, normal_cram, normal_crai, tumor_cram, tumor_crai, [], [], bed, tbi]
            }
        }

        RUN_STRELKA_SOMATIC(cram_pair_strelka,
                            fasta,
                            fasta_fai,
                            intervals_bed_combine_gz,
                            num_intervals)

        strelka_vcf = RUN_STRELKA_SOMATIC.out.strelka_vcf
        ch_versions = ch_versions.mix(RUN_STRELKA_SOMATIC.out.versions)
    }

    if (tools.contains('msisensorpro')) {
        MSISENSORPRO_MSI_SOMATIC(cram_pair_intervals, fasta, msisensorpro_scan)
        ch_versions = ch_versions.mix(MSISENSORPRO_MSI_SOMATIC.out.versions)
        msisensorpro_output = msisensorpro_output.mix(MSISENSORPRO_MSI_SOMATIC.out.output_report)
    }

    if (tools.contains('mutect2')) {
        // cram_pair_intervals.map{ meta, normal_cram, normal_crai, tumor_cram, tumor_crai, intervals ->
        //         [meta, [normal_cram, tumor_cram], [normal_crai, tumor_crai], intervals, ['normal']]
        //         }.set{cram_pair_mutect2}

        // GATK_TUMOR_NORMAL_SOMATIC_VARIANT_CALLING(
        //     cram_pair_mutect2,
        //     fasta,
        //     fasta_fai,
        //     dict,
        //     germline_resource,
        //     germline_resource_tbi,
        //     panel_of_normals,
        //     panel_of_normals_tbi,
        //     no_intervals,
        //     num_intervals,
        //     intervals_bed_combine_gz
        // )
        // ch_versions = ch_versions.mix(GATK_TUMOR_NORMAL_SOMATIC_VARIANT_CALLING.out.versions)
    }

    // if (tools.contains('tiddit')) {
    // }

    emit:
    manta_vcf
    msisensorpro_output
    mutect2_vcf
    strelka_vcf
    versions    = ch_versions
}
