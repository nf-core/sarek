//
// TUMOR VARIANT CALLING
// Should be only run on patients without normal sample
//

include { BAM_VARIANT_CALLING_CNVKIT                  } from '../bam_variant_calling_cnvkit/main'
include { BAM_VARIANT_CALLING_FREEBAYES               } from '../bam_variant_calling_freebayes/main'
include { BAM_VARIANT_CALLING_MPILEUP                 } from '../bam_variant_calling_mpileup/main'
include { BAM_VARIANT_CALLING_SINGLE_STRELKA          } from '../bam_variant_calling_single_strelka/main'
include { BAM_VARIANT_CALLING_SINGLE_TIDDIT           } from '../bam_variant_calling_single_tiddit/main'
include { BAM_VARIANT_CALLING_TUMOR_ONLY_CONTROLFREEC } from '../bam_variant_calling_tumor_only_controlfreec/main'
include { BAM_VARIANT_CALLING_TUMOR_ONLY_MANTA        } from '../bam_variant_calling_tumor_only_manta/main'
include { BAM_VARIANT_CALLING_TUMOR_ONLY_MUTECT2      } from '../bam_variant_calling_tumor_only_mutect2/main'

workflow BAM_VARIANT_CALLING_TUMOR_ONLY_ALL {
    take:
        tools                         // Mandatory, list of tools to apply
        cram_recalibrated             // channel: [mandatory] cram
        bwa                           // channel: [optional] bwa
        cf_chrom_len                  // channel: [optional] controlfreec length file
        chr_files
        cnvkit_reference
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
        panel_of_normals              // channel: [optional]  panel_of_normals
        panel_of_normals_tbi          // channel: [optional]  panel_of_normals_tbi

    main:

    ch_versions         = Channel.empty()

    //TODO: Temporary until the if's can be removed and printing to terminal is prevented with "when" in the modules.config
    freebayes_vcf       = Channel.empty()
    manta_vcf           = Channel.empty()
    mutect2_vcf         = Channel.empty()
    strelka_vcf         = Channel.empty()
    tiddit_vcf          = Channel.empty()

    // Remap channel with intervals
    cram_recalibrated_intervals = cram_recalibrated.combine(intervals)
        .map{ meta, cram, crai, intervals, num_intervals ->

            //If no interval file provided (0) then add empty list
            intervals_new = num_intervals == 0 ? [] : intervals

            [[
                data_type:      meta.data_type,
                id:             meta.sample,
                num_intervals:  num_intervals,
                patient:        meta.patient,
                sample:         meta.sample,
                sex:            meta.sex,
                status:         meta.status,
            ],
            cram, crai, intervals_new]
        }

    // Remap channel with gzipped intervals + indexes
    cram_recalibrated_intervals_gz_tbi = cram_recalibrated.combine(intervals_bed_gz_tbi)
        .map{ meta, cram, crai, bed_tbi, num_intervals ->

            //If no interval file provided (0) then add empty list
            bed_new = num_intervals == 0 ? [] : bed_tbi[0]
            tbi_new = num_intervals == 0 ? [] : bed_tbi[1]

            [[
                data_type:      meta.data_type,
                id:             meta.sample,
                num_intervals:  num_intervals,
                patient:        meta.patient,
                sample:         meta.sample,
                sex:            meta.sex,
                status:         meta.status,
            ],
            cram, crai, bed_new, tbi_new]
        }

    if (tools.split(',').contains('mpileup') || tools.split(',').contains('controlfreec')){
        cram_intervals_no_index = cram_recalibrated_intervals.map { meta, cram, crai, intervals ->
                                                                    [meta, cram, intervals]
                                                                    }
        BAM_VARIANT_CALLING_MPILEUP(
            cram_intervals_no_index,
            fasta
        )

        ch_versions = ch_versions.mix(BAM_VARIANT_CALLING_MPILEUP.out.versions)
    }

    if (tools.split(',').contains('controlfreec')){
        controlfreec_input = BAM_VARIANT_CALLING_MPILEUP.out.mpileup
                                .map{ meta, pileup_tumor ->
                                    [meta, [], pileup_tumor, [], [], [], []]
                                }

        length_file = cf_chrom_len ?: fasta_fai
        BAM_VARIANT_CALLING_TUMOR_ONLY_CONTROLFREEC(
            controlfreec_input,
            fasta,
            length_file,
            dbsnp,
            dbsnp_tbi,
            chr_files,
            mappability,
            intervals_bed_combined
        )

        ch_versions = ch_versions.mix(BAM_VARIANT_CALLING_TUMOR_ONLY_CONTROLFREEC.out.versions)
    }

    if(tools.split(',').contains('cnvkit')){
        cram_recalibrated_cnvkit_tumoronly = cram_recalibrated
            .map{ meta, cram, crai ->
                [meta, cram, []]
            }

        BAM_VARIANT_CALLING_CNVKIT (
            cram_recalibrated_cnvkit_tumoronly,
            fasta,
            fasta_fai,
            [],
            cnvkit_reference
        )

        ch_versions = ch_versions.mix(BAM_VARIANT_CALLING_CNVKIT.out.versions)
    }

    if (tools.split(',').contains('freebayes')){
        // Remap channel for Freebayes
        cram_recalibrated_intervals_freebayes = cram_recalibrated_intervals
            .map{ meta, cram, crai, intervals ->
                [meta, cram, crai, [], [], intervals]
            }

        BAM_VARIANT_CALLING_FREEBAYES(
            cram_recalibrated_intervals_freebayes,
            dict,
            fasta,
            fasta_fai
        )

        freebayes_vcf = BAM_VARIANT_CALLING_FREEBAYES.out.freebayes_vcf
        ch_versions   = ch_versions.mix(BAM_VARIANT_CALLING_FREEBAYES.out.versions)
    }

    if (tools.split(',').contains('mutect2')) {
        BAM_VARIANT_CALLING_TUMOR_ONLY_MUTECT2(
            cram_recalibrated_intervals,
            fasta,
            fasta_fai,
            dict,
            germline_resource,
            germline_resource_tbi,
            panel_of_normals,
            panel_of_normals_tbi
        )

        mutect2_vcf = BAM_VARIANT_CALLING_TUMOR_ONLY_MUTECT2.out.filtered_vcf
        ch_versions = ch_versions.mix(BAM_VARIANT_CALLING_TUMOR_ONLY_MUTECT2.out.versions)
    }

    if (tools.split(',').contains('manta')){

        BAM_VARIANT_CALLING_TUMOR_ONLY_MANTA(
            cram_recalibrated_intervals_gz_tbi,
            dict,
            fasta,
            fasta_fai
        )

        manta_vcf   = BAM_VARIANT_CALLING_TUMOR_ONLY_MANTA.out.manta_vcf
        ch_versions = ch_versions.mix(BAM_VARIANT_CALLING_TUMOR_ONLY_MANTA.out.versions)
    }

    if (tools.split(',').contains('strelka')) {

        BAM_VARIANT_CALLING_SINGLE_STRELKA(
            cram_recalibrated_intervals_gz_tbi,
            dict,
            fasta,
            fasta_fai
        )

        strelka_vcf = BAM_VARIANT_CALLING_SINGLE_STRELKA.out.strelka_vcf
        ch_versions = ch_versions.mix(BAM_VARIANT_CALLING_SINGLE_STRELKA.out.versions)
    }

        //TIDDIT
    if (tools.split(',').contains('tiddit')){

        BAM_VARIANT_CALLING_SINGLE_TIDDIT(
            cram_recalibrated,
            fasta.map{ it -> [[id:it[0].baseName], it] },
            bwa
        )

        tiddit_vcf = BAM_VARIANT_CALLING_SINGLE_TIDDIT.out.tiddit_vcf
        ch_versions = ch_versions.mix(BAM_VARIANT_CALLING_SINGLE_TIDDIT.out.versions)
    }


    emit:
    freebayes_vcf
    manta_vcf
    mutect2_vcf
    strelka_vcf
    tiddit_vcf

    versions = ch_versions
}
