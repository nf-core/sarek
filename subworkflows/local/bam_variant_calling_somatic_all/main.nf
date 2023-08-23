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
include { MSISENSORPRO_MSISOMATIC                       } from '../../../modules/nf-core/msisensorpro/msisomatic/main'

workflow BAM_VARIANT_CALLING_SOMATIC_ALL {
    take:
    tools                         // Mandatory, list of tools to apply
    cram                          // channel: [mandatory] cram
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
    intervals                     // channel: [mandatory] [ intervals, num_intervals ] or [ [], 0 ] if no intervals
    intervals_bed_gz_tbi          // channel: [mandatory] intervals/target regions index zipped and indexed
    intervals_bed_combined        // channel: [mandatory] intervals/target regions in one file unzipped
    intervals_bed_gz_tbi_combined // channel: [mandatory] intervals/target regions in one file zipped
    mappability
    msisensorpro_scan             // channel: [optional]  msisensorpro_scan
    panel_of_normals              // channel: [optional]  panel_of_normals
    panel_of_normals_tbi          // channel: [optional]  panel_of_normals_tbi
    allele_files                  // channel: [optional]  ascat allele files
    loci_files                    // channel: [optional]  ascat loci files
    gc_file                       // channel: [optional]  ascat gc content file
    rt_file                       // channel: [optional]  ascat rt file
    joint_mutect2                 // boolean: [mandatory] [default: false] run mutect2 in joint mode

    main:
    versions          = Channel.empty()

    //TODO: Temporary until the if's can be removed and printing to terminal is prevented with "when" in the modules.config
    vcf_freebayes       = Channel.empty()
    vcf_manta           = Channel.empty()
    vcf_strelka         = Channel.empty()
    out_msisensorpro    = Channel.empty()
    vcf_mutect2         = Channel.empty()
    vcf_tiddit          = Channel.empty()

    if (checkTools(params.tools, 'ascat')) {
        BAM_VARIANT_CALLING_SOMATIC_ASCAT(
            cram,
            allele_files,
            loci_files,
            intervals_bed_combined,
            fasta,
            gc_file,
            rt_file
        )

        versions = versions.mix(BAM_VARIANT_CALLING_SOMATIC_ASCAT.out.versions)
    }

    // CONTROLFREEC
    if (checkTools(params.tools, 'controlfreec')) {
        // Remap channels to match module/subworkflow
        cram_normal = cram.map { meta, normal_cram, normal_crai, tumor_cram, tumor_crai -> [ meta, normal_cram, normal_crai ] }
        cram_tumor = cram.map { meta, normal_cram, normal_crai, tumor_cram, tumor_crai -> [ meta, tumor_cram, tumor_crai ] }

        MPILEUP_NORMAL(
            cram_normal,
            dict,
            fasta,
            intervals
        )

        MPILEUP_TUMOR(
            cram_tumor,
            dict,
            fasta,
            intervals
        )

        mpileup_normal = MPILEUP_NORMAL.out.mpileup
        mpileup_tumor = MPILEUP_TUMOR.out.mpileup
            // Remap channel to match module/subworkflow
        mpileup_pair = mpileup_normal.cross(mpileup_tumor).map{ normal, tumor -> [ normal[0], normal[1], tumor[1], [], [], [], [] ] }

        length_file = cf_chrom_len ?: fasta_fai

        BAM_VARIANT_CALLING_SOMATIC_CONTROLFREEC(
            mpileup_pair,
            fasta,
            length_file,
            dbsnp,
            dbsnp_tbi,
            chr_files,
            mappability,
            intervals_bed_combined
        )

        versions = versions.mix(MPILEUP_NORMAL.out.versions)
        versions = versions.mix(MPILEUP_TUMOR.out.versions)
        versions = versions.mix(BAM_VARIANT_CALLING_SOMATIC_CONTROLFREEC.out.versions)
    }

    // CNVKIT
    if (checkTools(params.tools, 'cnvkit')) {
        BAM_VARIANT_CALLING_CNVKIT(
            // Remap channel to match module/subworkflow
            cram.map{ meta, normal_cram, normal_crai, tumor_cram, tumor_crai -> [ meta, tumor_cram, normal_cram ] },
            fasta,
            fasta_fai,
            intervals_bed_combined,
            []
        )

        versions = versions.mix(BAM_VARIANT_CALLING_CNVKIT.out.versions)
    }

    // FREEBAYES
    if (checkTools(params.tools, 'freebayes')) {
        BAM_VARIANT_CALLING_FREEBAYES(
            cram,
            dict,
            fasta,
            fasta_fai,
            intervals
        )

        vcf_freebayes = BAM_VARIANT_CALLING_FREEBAYES.out.vcf
        versions   = versions.mix(BAM_VARIANT_CALLING_FREEBAYES.out.versions)
    }

    // MANTA
    if (checkTools(params.tools, 'manta')) {
        BAM_VARIANT_CALLING_SOMATIC_MANTA(
            cram,
            fasta,
            fasta_fai,
            intervals_bed_gz_tbi_combined
        )

        vcf_manta = BAM_VARIANT_CALLING_SOMATIC_MANTA.out.vcf
        versions = versions.mix(BAM_VARIANT_CALLING_SOMATIC_MANTA.out.versions)
    }

    // STRELKA
    if (checkTools(params.tools, 'strelka')) {
        // Remap channel to match module/subworkflow
        cram_strelka = (checkTools(params.tools, 'manta')) ?
            cram.join(BAM_VARIANT_CALLING_SOMATIC_MANTA.out.candidate_small_indels_vcf, failOnDuplicate: true, failOnMismatch: true).join(BAM_VARIANT_CALLING_SOMATIC_MANTA.out.candidate_small_indels_vcf_tbi, failOnDuplicate: true, failOnMismatch: true) :
            cram.map{ meta, normal_cram, normal_crai, tumor_cram, tumor_crai -> [ meta, normal_cram, normal_crai, tumor_cram, tumor_crai, [], [] ] }

        BAM_VARIANT_CALLING_SOMATIC_STRELKA(
            cram_strelka,
            // Remap channel to match module/subworkflow
            dict,
            fasta,
            fasta_fai,
            intervals_bed_gz_tbi
        )

        vcf_strelka = Channel.empty().mix(BAM_VARIANT_CALLING_SOMATIC_STRELKA.out.vcf)
        versions = versions.mix(BAM_VARIANT_CALLING_SOMATIC_STRELKA.out.versions)
    }

    // MSISENSOR
    if (checkTools(params.tools, 'msisensorpro')) {
        MSISENSORPRO_MSISOMATIC(cram.combine(intervals_bed_combined), fasta, msisensorpro_scan)

        versions = versions.mix(MSISENSORPRO_MSISOMATIC.out.versions)
        out_msisensorpro = out_msisensorpro.mix(MSISENSORPRO_MSISOMATIC.out.output_report)
    }

    // MUTECT2
    if (checkTools(params.tools, 'mutect2')) {
        BAM_VARIANT_CALLING_SOMATIC_MUTECT2(
            // Remap channel to match module/subworkflow
            cram.map { meta, normal_cram, normal_crai, tumor_cram, tumor_crai -> [ meta, [ normal_cram, tumor_cram ], [ normal_crai, tumor_crai ] ] },
            // Remap channel to match module/subworkflow
            fasta.map{ it -> [ [ id:'fasta' ], it ] },
            // Remap channel to match module/subworkflow
            fasta_fai.map{ it -> [ [ id:'fasta_fai' ], it ] },
            dict,
            germline_resource,
            germline_resource_tbi,
            panel_of_normals,
            panel_of_normals_tbi,
            intervals,
            joint_mutect2
        )

        vcf_mutect2 = BAM_VARIANT_CALLING_SOMATIC_MUTECT2.out.vcf_filtered
        versions = versions.mix(BAM_VARIANT_CALLING_SOMATIC_MUTECT2.out.versions)
    }

    // TIDDIT
    if (checkTools(params.tools, 'tiddit')) {
        BAM_VARIANT_CALLING_SOMATIC_TIDDIT(
            // Remap channel to match module/subworkflow
            cram.map{ meta, normal_cram, normal_crai, tumor_cram, tumor_crai -> [ meta, normal_cram, normal_crai ] },
            // Remap channel to match module/subworkflow
            cram.map{ meta, normal_cram, normal_crai, tumor_cram, tumor_crai -> [ meta, tumor_cram, tumor_crai ] },
            // Remap channel to match module/subworkflow
            fasta.map{ it -> [ [ id:'fasta' ], it ] },
            bwa)
        vcf_tiddit = BAM_VARIANT_CALLING_SOMATIC_TIDDIT.out.vcf
        versions = versions.mix(BAM_VARIANT_CALLING_SOMATIC_TIDDIT.out.versions)
    }

    vcf_all = Channel.empty().mix(
        vcf_freebayes,
        vcf_manta,
        vcf_mutect2,
        vcf_strelka,
        vcf_tiddit
    )

    emit:
    out_msisensorpro
    vcf_all
    vcf_freebayes
    vcf_manta
    vcf_mutect2
    vcf_strelka
    vcf_tiddit

    versions
}
