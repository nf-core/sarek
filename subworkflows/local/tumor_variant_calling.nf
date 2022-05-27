//
// TUMOR VARIANT CALLING
// Should be only run on patients without normal sample
//

include { RUN_FREEBAYES                           } from '../nf-core/variantcalling/freebayes/main.nf'
include { GATK_TUMOR_ONLY_SOMATIC_VARIANT_CALLING } from '../../subworkflows/nf-core/gatk4/tumor_only_somatic_variant_calling/main'
include { RUN_MANTA_TUMORONLY                     } from '../nf-core/variantcalling/manta/tumoronly/main.nf'
include { RUN_STRELKA_SINGLE                      } from '../nf-core/variantcalling/strelka/single/main.nf'
include { RUN_CONTROLFREEC_TUMORONLY              } from '../nf-core/variantcalling/controlfreec/tumoronly/main.nf'

workflow TUMOR_ONLY_VARIANT_CALLING {
    take:
        tools                        // Mandatory, list of tools to apply
        cram_recalibrated            // channel: [mandatory] cram
        dbsnp                        // channel: [mandatory] dbsnp
        dbsnp_tbi                    // channel: [mandatory] dbsnp_tbi
        dict                         // channel: [mandatory] dict
        fasta                        // channel: [mandatory] fasta
        fasta_fai                    // channel: [mandatory] fasta_fai
        intervals                    // channel: [mandatory] intervals/target regions
        intervals_bed_gz_tbi         // channel: [mandatory] intervals/target regions index zipped and indexed
        intervals_bed_combine_gz_tbi // channel: [mandatory] intervals/target regions index zipped and indexed
        intervals_bed_combine_gz     // channel: [mandatory] intervals/target regions index zipped and indexed in one file
        intervals_bed_combined        // channel: [mandatory] intervals/target regions in one file unzipped
        germline_resource            // channel: [optional]  germline_resource
        germline_resource_tbi        // channel: [optional]  germline_resource_tbi
        panel_of_normals             // channel: [optional]  panel_of_normals
        panel_of_normals_tbi         // channel: [optional]  panel_of_normals_tbi
        chr_files
        mappability

    main:

    ch_versions         = Channel.empty()

    //TODO: Temporary until the if's can be removed and printing to terminal is prevented with "when" in the modules.config
    freebayes_vcf       = Channel.empty()
    manta_vcf           = Channel.empty()
    mutect2_vcf         = Channel.empty()
    strelka_vcf         = Channel.empty()

    // Remap channel with intervals
    cram_recalibrated_intervals = cram_recalibrated.combine(intervals)
        .map{ meta, cram, crai, intervals, num_intervals ->

            //If no interval file provided (0) then add empty list
            intervals_new = num_intervals == 0 ? [] : intervals

            // If either no scatter/gather is done, i.e. no interval (0) or one interval (1), then don't rename samples
            new_id = num_intervals <= 1 ? meta.sample : meta.sample + "_" + intervals_new.baseName

            [[patient:meta.patient, sample:meta.sample, gender:meta.gender, status:meta.status, id:new_id, data_type:meta.data_type, num_intervals:num_intervals],
            cram, crai, intervals_new]
        }

    // Remap channel with gzipped intervals + indexes
    cram_recalibrated_intervals_gz_tbi = cram_recalibrated.combine(intervals_bed_gz_tbi)
        .map{ meta, cram, crai, bed_tbi, num_intervals ->

            //If no interval file provided (0) then add empty list
            bed_new = num_intervals == 0 ? [] : bed_tbi[0]
            tbi_new = num_intervals == 0 ? [] : bed_tbi[1]

            // If either no scatter/gather is done, i.e. no interval (0) or one interval (1), then don't rename samples
            new_id = num_intervals <= 1 ? meta.sample : meta.sample + "_" + bed_new.simpleName

            [[patient:meta.patient, sample:meta.sample, gender:meta.gender, status:meta.status, id:new_id, data_type:meta.data_type, num_intervals:num_intervals],
            cram, crai, bed_new, tbi_new]
        }

    if(tools.contains('controlfreec')){
        cram_intervals_no_index = cram_recalibrated_intervals.map { meta, cram, crai, intervals ->
                                                                    [meta, cram, intervals]
                                                                    }
        RUN_CONTROLFREEC_TUMORONLY(
                        cram_intervals_no_index,
                        fasta,
                        fasta_fai,
                        dbsnp,
                        dbsnp_tbi,
                        chr_files,
                        mappability,
                        intervals_bed_combined)
        ch_versions = ch_versions.mix(RUN_CONTROLFREEC_TUMORONLY.out.versions)
    }

    if (tools.contains('freebayes')){
        // Remap channel for Freebayes
        cram_recalibrated_intervals_freebayes = cram_recalibrated_intervals
            .map{ meta, cram, crai, intervals ->
                [meta, cram, crai, [], [], intervals]
            }

        RUN_FREEBAYES(cram_recalibrated_intervals_freebayes, fasta, fasta_fai, intervals_bed_combine_gz)

        freebayes_vcf = RUN_FREEBAYES.out.freebayes_vcf
        ch_versions   = ch_versions.mix(RUN_FREEBAYES.out.versions)
    }

    if (tools.contains('mutect2')) {
        GATK_TUMOR_ONLY_SOMATIC_VARIANT_CALLING(cram_recalibrated_intervals,
                                                fasta,
                                                fasta_fai,
                                                dict,
                                                germline_resource,
                                                germline_resource_tbi,
                                                panel_of_normals,
                                                panel_of_normals_tbi,
                                                intervals_bed_combine_gz)

        mutect2_vcf = GATK_TUMOR_ONLY_SOMATIC_VARIANT_CALLING.out.filtered_vcf
        ch_versions = ch_versions.mix(GATK_TUMOR_ONLY_SOMATIC_VARIANT_CALLING.out.versions)
    }

    if (tools.contains('manta')){
        RUN_MANTA_TUMORONLY(cram_recalibrated_intervals_gz_tbi,
                            fasta,
                            fasta_fai,
                            intervals_bed_combine_gz)

        manta_vcf   = RUN_MANTA_TUMORONLY.out.manta_vcf
        ch_versions = ch_versions.mix(RUN_MANTA_TUMORONLY.out.versions)
    }

    if (tools.contains('strelka')) {
        RUN_STRELKA_SINGLE( cram_recalibrated_intervals_gz_tbi,
                            fasta,
                            fasta_fai,
                            intervals_bed_combine_gz)

        strelka_vcf = RUN_STRELKA_SINGLE.out.strelka_vcf
        ch_versions = ch_versions.mix(RUN_STRELKA_SINGLE.out.versions)
    }


    // if (tools.contains('tiddit')){
    // }

    emit:
    versions = ch_versions

    freebayes_vcf
    manta_vcf
    mutect2_vcf
    strelka_vcf

}
