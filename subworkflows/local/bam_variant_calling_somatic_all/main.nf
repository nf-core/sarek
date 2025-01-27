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
include { BAM_VARIANT_CALLING_INDEXCOV                  } from '../bam_variant_calling_indexcov/main'
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
    wes                           // boolean: [mandatory] [default: false] whether targeted data is processed

    main:
    versions          = Channel.empty()

    //TODO: Temporary until the if's can be removed and printing to terminal is prevented with "when" in the modules.config
    vcf_freebayes       = Channel.empty()
    vcf_manta           = Channel.empty()
    vcf_strelka         = Channel.empty()
    out_msisensorpro    = Channel.empty()
    vcf_mutect2         = Channel.empty()
    vcf_tiddit          = Channel.empty()
    out_indexcov        = Channel.empty()

    if (tools.split(',').contains('ascat')) {
        BAM_VARIANT_CALLING_SOMATIC_ASCAT(
            cram,
            allele_files,
            loci_files,
            (wes ? intervals_bed_combined : []), // No intervals needed if not WES
            fasta.map{ meta, fasta -> [ fasta ] },
            gc_file,
            rt_file
        )

        versions = versions.mix(BAM_VARIANT_CALLING_SOMATIC_ASCAT.out.versions)
    }

    // CONTROLFREEC
    if (tools.split(',').contains('controlfreec')) {
        // Remap channels to match module/subworkflow
        cram_normal = cram.map{ meta, normal_cram, normal_crai, tumor_cram, tumor_crai -> [ meta, normal_cram, normal_crai ] }
        cram_tumor = cram.map{ meta, normal_cram, normal_crai, tumor_cram, tumor_crai -> [ meta, tumor_cram, tumor_crai ] }

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

        length_file  = cf_chrom_len ?: fasta_fai.map{ meta, fasta_fai -> [ fasta_fai ] }

        intervals_controlfreec = wes ? intervals_bed_combined : []

        BAM_VARIANT_CALLING_SOMATIC_CONTROLFREEC(
            mpileup_pair,
            fasta.map{ meta, fasta -> [ fasta ] },
            length_file,
            dbsnp,
            dbsnp_tbi,
            chr_files,
            mappability,
            intervals_controlfreec
        )

        versions = versions.mix(MPILEUP_NORMAL.out.versions)
        versions = versions.mix(MPILEUP_TUMOR.out.versions)
        versions = versions.mix(BAM_VARIANT_CALLING_SOMATIC_CONTROLFREEC.out.versions)
    }

    // CNVKIT
    if (tools.split(',').contains('cnvkit')) {
        BAM_VARIANT_CALLING_CNVKIT(
            // Remap channel to match module/subworkflow
            cram.map{ meta, normal_cram, normal_crai, tumor_cram, tumor_crai -> [ meta, tumor_cram, normal_cram ] },
            fasta,
            fasta_fai,
            intervals_bed_combined.map{ it -> [[id:it[0].baseName], it] },
            [[id:"null"], []]
        )

        versions = versions.mix(BAM_VARIANT_CALLING_CNVKIT.out.versions)
    }

    // FREEBAYES
    if (tools.split(',').contains('freebayes')) {
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
    if (tools.split(',').contains('manta')) {
        BAM_VARIANT_CALLING_SOMATIC_MANTA(
            cram,
            fasta,
            fasta_fai,
            intervals_bed_gz_tbi_combined
        )

        vcf_manta = BAM_VARIANT_CALLING_SOMATIC_MANTA.out.vcf
        versions = versions.mix(BAM_VARIANT_CALLING_SOMATIC_MANTA.out.versions)
    }


    // INDEXCOV, for WGS only
    if (params.wes==false &&  tools.split(',').contains('indexcov')) {
        BAM_VARIANT_CALLING_INDEXCOV (
            cram,
            fasta,
            fasta_fai
        )

        out_indexcov = BAM_VARIANT_CALLING_INDEXCOV.out.out_indexcov
        versions = versions.mix(BAM_VARIANT_CALLING_INDEXCOV.out.versions)
    }


    // STRELKA
    if (tools.split(',').contains('strelka')) {
        // Remap channel to match module/subworkflow
        cram_strelka = (tools.split(',').contains('manta')) ?
            cram.join(BAM_VARIANT_CALLING_SOMATIC_MANTA.out.candidate_small_indels_vcf, failOnDuplicate: true, failOnMismatch: true).join(BAM_VARIANT_CALLING_SOMATIC_MANTA.out.candidate_small_indels_vcf_tbi, failOnDuplicate: true, failOnMismatch: true) :
            cram.map{ meta, normal_cram, normal_crai, tumor_cram, tumor_crai -> [ meta, normal_cram, normal_crai, tumor_cram, tumor_crai, [], [] ] }

        BAM_VARIANT_CALLING_SOMATIC_STRELKA(
            cram_strelka,
            // Remap channel to match module/subworkflow
            dict,
            fasta.map{ meta, fasta -> [ fasta ] },
            fasta_fai.map{ meta, fasta_fai -> [ fasta_fai ] },
            intervals_bed_gz_tbi
        )

        vcf_strelka = Channel.empty().mix(BAM_VARIANT_CALLING_SOMATIC_STRELKA.out.vcf)
        versions = versions.mix(BAM_VARIANT_CALLING_SOMATIC_STRELKA.out.versions)
    }

    // MSISENSOR
    if (tools.split(',').contains('msisensorpro')) {
        MSISENSORPRO_MSISOMATIC(cram.combine(intervals_bed_combined), fasta.map{ meta, fasta -> [ fasta ] }, msisensorpro_scan)

        versions = versions.mix(MSISENSORPRO_MSISOMATIC.out.versions)
        out_msisensorpro = out_msisensorpro.mix(MSISENSORPRO_MSISOMATIC.out.output_report)
    }

    // MUTECT2
    if (tools.split(',').contains('mutect2')) {
        BAM_VARIANT_CALLING_SOMATIC_MUTECT2(
            // Remap channel to match module/subworkflow
            // Adjust meta.map to simplify joining channels
            // joint_mutect2 mode needs different meta.map than regular mode
            cram.map{ meta, normal_cram, normal_crai, tumor_cram, tumor_crai ->
                joint_mutect2 ?
                //we need to keep all fields and then remove on a per-tool-basis to ensure proper joining at the filtering step
                [ meta + [ id:meta.patient ], [ normal_cram, tumor_cram ], [ normal_crai, tumor_crai ] ] :
                [ meta, [ normal_cram, tumor_cram ], [ normal_crai, tumor_crai ] ]
            },
            fasta,
            fasta_fai,
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
    if (tools.split(',').contains('tiddit')) {
        BAM_VARIANT_CALLING_SOMATIC_TIDDIT(
            // Remap channel to match module/subworkflow
            cram.map{ meta, normal_cram, normal_crai, tumor_cram, tumor_crai -> [ meta, normal_cram, normal_crai ] },
            // Remap channel to match module/subworkflow
            cram.map{ meta, normal_cram, normal_crai, tumor_cram, tumor_crai -> [ meta, tumor_cram, tumor_crai ] },
            fasta,
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
    out_indexcov
    out_msisensorpro
    vcf_all
    vcf_freebayes
    vcf_manta
    vcf_mutect2
    vcf_strelka
    vcf_tiddit

    versions
}
