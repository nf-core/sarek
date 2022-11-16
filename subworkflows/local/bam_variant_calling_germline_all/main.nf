//
// GERMLINE VARIANT CALLING
//

include { BAM_VARIANT_CALLING_CNVKIT          } from '../bam_variant_calling_cnvkit/main'
include { BAM_VARIANT_CALLING_DEEPVARIANT     } from '../bam_variant_calling_deepvariant/main'
include { BAM_VARIANT_CALLING_FREEBAYES       } from '../bam_variant_calling_freebayes/main'
include { BAM_VARIANT_CALLING_GERMLINE_MANTA  } from '../bam_variant_calling_germline_manta/main'
include { BAM_VARIANT_CALLING_HAPLOTYPECALLER } from '../bam_variant_calling_haplotypecaller/main'
include { BAM_VARIANT_CALLING_MPILEUP         } from '../bam_variant_calling_mpileup/main'
include { BAM_VARIANT_CALLING_SINGLE_STRELKA  } from '../bam_variant_calling_single_strelka/main'
include { BAM_VARIANT_CALLING_SINGLE_TIDDIT   } from '../bam_variant_calling_single_tiddit/main'

workflow BAM_VARIANT_CALLING_GERMLINE_ALL {
    take:
        tools                             // Mandatory, list of tools to apply
        cram_recalibrated                 // channel: [mandatory] cram
        bwa                               // channel: [mandatory] bwa
        dbsnp                             // channel: [mandatory] dbsnp
        dbsnp_tbi                         // channel: [mandatory] dbsnp_tbi
        dict                              // channel: [mandatory] dict
        fasta                             // channel: [mandatory] fasta
        fasta_fai                         // channel: [mandatory] fasta_fai
        intervals                         // channel: [mandatory] intervals/target regions
        intervals_bed_gz_tbi              // channel: [mandatory] intervals/target regions index zipped and indexed
        intervals_bed_combined            // channel: [mandatory] intervals/target regions in one file unzipped
        intervals_bed_combined_haplotypec // channel: [mandatory] intervals/target regions in one file unzipped, no_intervals.bed if no_intervals
        known_sites_indels
        known_sites_indels_tbi
        known_sites_snps
        known_sites_snps_tbi

    main:

    ch_versions         = Channel.empty()

    //TODO: Temporary until the if's can be removed and printing to terminal is prevented with "when" in the modules.config
    deepvariant_vcf     = Channel.empty()
    freebayes_vcf       = Channel.empty()
    genotype_gvcf       = Channel.empty()
    haplotypecaller_vcf = Channel.empty()
    manta_vcf           = Channel.empty()
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

    if(tools.split(',').contains('mpileup')){
        cram_intervals_no_index = cram_recalibrated_intervals
            .map { meta, cram, crai, intervals ->
                [meta, cram, intervals]
            }

        BAM_VARIANT_CALLING_MPILEUP(
            cram_intervals_no_index,
            fasta
        )

        mpileup_germline = BAM_VARIANT_CALLING_MPILEUP.out.mpileup
        ch_versions = ch_versions.mix(BAM_VARIANT_CALLING_MPILEUP.out.versions)
    }

    // CNVKIT

    if(tools.split(',').contains('cnvkit')){
        cram_recalibrated_cnvkit_germline = cram_recalibrated
            .map{ meta, cram, crai ->
                [meta, [], cram]
            }

        BAM_VARIANT_CALLING_CNVKIT(
            cram_recalibrated_cnvkit_germline,
            fasta,
            fasta_fai,
            intervals_bed_combined,
            []
        )

        ch_versions     = ch_versions.mix(BAM_VARIANT_CALLING_CNVKIT.out.versions)
    }

    // DEEPVARIANT
    if(tools.split(',').contains('deepvariant')){
        BAM_VARIANT_CALLING_DEEPVARIANT(
            cram_recalibrated_intervals,
            dict,
            fasta,
            fasta_fai
        )

        deepvariant_vcf = Channel.empty().mix(BAM_VARIANT_CALLING_DEEPVARIANT.out.deepvariant_vcf)
        ch_versions     = ch_versions.mix(BAM_VARIANT_CALLING_DEEPVARIANT.out.versions)
    }

    // FREEBAYES
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

        freebayes_vcf   = BAM_VARIANT_CALLING_FREEBAYES.out.freebayes_vcf
        ch_versions     = ch_versions.mix(BAM_VARIANT_CALLING_FREEBAYES.out.versions)
    }

    // HAPLOTYPECALLER
    if (tools.split(',').contains('haplotypecaller')){
        cram_recalibrated_intervals_haplotypecaller = cram_recalibrated_intervals
            .map{ meta, cram, crai, intervals ->

            intervals_name = meta.num_intervals == 0 ? "no_interval" : intervals.simpleName
            new_meta = params.joint_germline ? [
                                                    data_type:meta.data_type,
                                                    id:meta.sample,
                                                    intervals_name:intervals_name,
                                                    num_intervals:meta.num_intervals,
                                                    patient:meta.patient,
                                                    sample:meta.sample,
                                                    sex:meta.sex,
                                                    status:meta.status
                                                ]
                                            : meta

                [new_meta, cram, crai, intervals, []]
        }

        BAM_VARIANT_CALLING_HAPLOTYPECALLER(cram_recalibrated_intervals_haplotypecaller,
                        fasta,
                        fasta_fai,
                        dict,
                        dbsnp,
                        dbsnp_tbi,
                        known_sites_indels,
                        known_sites_indels_tbi,
                        known_sites_snps,
                        known_sites_snps_tbi,
                        intervals_bed_combined_haplotypec)

        haplotypecaller_vcf  = BAM_VARIANT_CALLING_HAPLOTYPECALLER.out.filtered_vcf
        ch_versions          = ch_versions.mix(BAM_VARIANT_CALLING_HAPLOTYPECALLER.out.versions)
    }

    // MANTA
    if (tools.split(',').contains('manta')){
        BAM_VARIANT_CALLING_GERMLINE_MANTA (
            cram_recalibrated_intervals_gz_tbi,
            dict,
            fasta,
            fasta_fai
        )

        manta_vcf   = BAM_VARIANT_CALLING_GERMLINE_MANTA.out.manta_vcf
        ch_versions = ch_versions.mix(BAM_VARIANT_CALLING_GERMLINE_MANTA.out.versions)
    }

    // STRELKA
    if (tools.split(',').contains('strelka')){
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
    deepvariant_vcf
    freebayes_vcf
    genotype_gvcf
    haplotypecaller_vcf
    manta_vcf
    strelka_vcf
    tiddit_vcf

    versions = ch_versions
}
