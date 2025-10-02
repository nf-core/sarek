//
// PAIRED VARIANT CALLING
//

include { BAM_VARIANT_CALLING_CNVKIT                    } from '../bam_variant_calling_cnvkit'
include { BAM_VARIANT_CALLING_FREEBAYES                 } from '../bam_variant_calling_freebayes'
include { BAM_VARIANT_CALLING_INDEXCOV                  } from '../bam_variant_calling_indexcov'
include { BAM_VARIANT_CALLING_MPILEUP as MPILEUP_NORMAL } from '../bam_variant_calling_mpileup'
include { BAM_VARIANT_CALLING_MPILEUP as MPILEUP_TUMOR  } from '../bam_variant_calling_mpileup'
include { BAM_VARIANT_CALLING_SOMATIC_ASCAT             } from '../bam_variant_calling_somatic_ascat'
include { BAM_VARIANT_CALLING_SOMATIC_CONTROLFREEC      } from '../bam_variant_calling_somatic_controlfreec'
include { BAM_VARIANT_CALLING_SOMATIC_MANTA             } from '../bam_variant_calling_somatic_manta'
include { BAM_VARIANT_CALLING_SOMATIC_MUSE              } from '../bam_variant_calling_somatic_muse'
include { BAM_VARIANT_CALLING_SOMATIC_MUTECT2           } from '../bam_variant_calling_somatic_mutect2'
include { BAM_VARIANT_CALLING_SOMATIC_STRELKA           } from '../bam_variant_calling_somatic_strelka'
include { BAM_VARIANT_CALLING_SOMATIC_TIDDIT            } from '../bam_variant_calling_somatic_tiddit'
include { BAM_VARIANT_CALLING_SOMATIC_TNSCOPE           } from '../bam_variant_calling_somatic_tnscope'
include { MSISENSOR2_MSI                                } from '../../../modules/nf-core/msisensor2/msi'
include { MSISENSORPRO_MSISOMATIC                       } from '../../../modules/nf-core/msisensorpro/msisomatic'

workflow BAM_VARIANT_CALLING_SOMATIC_ALL {
    take:
    tools                         // Mandatory, list of tools to apply
    bam                           // channel: [mandatory] bam
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
    msisensor2_scan               // channel: [optional]  msisensor2_scan
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
    // Channels are often remapped to match module/subworkflow

    // Gather all versions
    versions = Channel.empty()

    //TODO: Temporary until the if's can be removed and printing to terminal is prevented with "when" in the modules.config
    out_indexcov     = Channel.empty()
    out_msisensor2   = Channel.empty()
    out_msisensorpro = Channel.empty()
    vcf_freebayes    = Channel.empty()
    vcf_manta        = Channel.empty()
    vcf_muse         = Channel.empty()
    vcf_mutect2      = Channel.empty()
    vcf_strelka      = Channel.empty()
    vcf_tiddit       = Channel.empty()
    vcf_tnscope      = Channel.empty()

    if (tools && tools.split(',').contains('ascat')) {
        BAM_VARIANT_CALLING_SOMATIC_ASCAT(
            cram,
            allele_files,
            loci_files,
            (wes ? intervals_bed_combined : []),
            fasta.map { _meta, fasta_ -> [fasta_] },
            gc_file,
            rt_file,
        )

        versions = versions.mix(BAM_VARIANT_CALLING_SOMATIC_ASCAT.out.versions)
    }

    // CONTROLFREEC
    if (tools && tools.split(',').contains('controlfreec')) {
        cram_normal = cram.map { meta, normal_cram, normal_crai, _tumor_cram, _tumor_crai -> [meta, normal_cram, normal_crai] }
        cram_tumor = cram.map { meta, _normal_cram, _normal_crai, tumor_cram, tumor_crai -> [meta, tumor_cram, tumor_crai] }

        MPILEUP_NORMAL(
            cram_normal,
            dict,
            fasta,
            intervals,
        )

        MPILEUP_TUMOR(
            cram_tumor,
            dict,
            fasta,
            intervals,
        )

        mpileup_normal = MPILEUP_NORMAL.out.mpileup
        mpileup_tumor = MPILEUP_TUMOR.out.mpileup
        mpileup_pair = mpileup_normal.cross(mpileup_tumor).map { normal, tumor -> [normal[0], normal[1], tumor[1], [], [], [], []] }

        BAM_VARIANT_CALLING_SOMATIC_CONTROLFREEC(
            mpileup_pair,
            fasta.map { _meta, fasta_ -> [fasta_] },
            cf_chrom_len ?: fasta_fai.map { _meta, fasta_fai_ -> [fasta_fai_] },
            dbsnp,
            dbsnp_tbi,
            chr_files,
            mappability,
            wes ? intervals_bed_combined : [],
        )

        versions = versions.mix(MPILEUP_NORMAL.out.versions)
        versions = versions.mix(MPILEUP_TUMOR.out.versions)
        versions = versions.mix(BAM_VARIANT_CALLING_SOMATIC_CONTROLFREEC.out.versions)
    }

    // CNVKIT
    if (tools && tools.split(',').contains('cnvkit')) {
        BAM_VARIANT_CALLING_CNVKIT(
            bam.map { meta, normal_bam, _normal_bai, tumor_bam, _tumor_bai -> [meta, tumor_bam, normal_bam] },
            fasta,
            fasta_fai,
            intervals_bed_combined.map { _intervals -> _intervals ? [[id: _intervals[0].baseName], _intervals] : [[id: 'no_intervals'], []] },
            [[id: "null"], []],
        )

        versions = versions.mix(BAM_VARIANT_CALLING_CNVKIT.out.versions)
    }

    // FREEBAYES
    if (tools && tools.split(',').contains('freebayes')) {
        BAM_VARIANT_CALLING_FREEBAYES(
            cram,
            dict,
            fasta,
            fasta_fai,
            intervals,
        )

        vcf_freebayes = BAM_VARIANT_CALLING_FREEBAYES.out.vcf
        versions = versions.mix(BAM_VARIANT_CALLING_FREEBAYES.out.versions)
    }

    // MANTA
    if (tools && tools.split(',').contains('manta')) {
        BAM_VARIANT_CALLING_SOMATIC_MANTA(
            cram,
            fasta,
            fasta_fai,
            intervals_bed_gz_tbi_combined,
        )

        vcf_manta = BAM_VARIANT_CALLING_SOMATIC_MANTA.out.vcf
        versions = versions.mix(BAM_VARIANT_CALLING_SOMATIC_MANTA.out.versions)
    }


    // INDEXCOV
    //   WGS only
    if (params.wes == false && tools.split(',').contains('indexcov')) {
        BAM_VARIANT_CALLING_INDEXCOV(
            cram,
            fasta,
            fasta_fai,
        )

        out_indexcov = BAM_VARIANT_CALLING_INDEXCOV.out.out_indexcov
        versions = versions.mix(BAM_VARIANT_CALLING_INDEXCOV.out.versions)
    }


    // STRELKA
    if (tools && tools.split(',').contains('strelka')) {
        cram_strelka = tools.split(',').contains('manta')
            ? cram.join(BAM_VARIANT_CALLING_SOMATIC_MANTA.out.candidate_small_indels_vcf, failOnDuplicate: true, failOnMismatch: true).join(BAM_VARIANT_CALLING_SOMATIC_MANTA.out.candidate_small_indels_vcf_tbi, failOnDuplicate: true, failOnMismatch: true)
            : cram.map { meta, normal_cram, normal_crai, tumor_cram, tumor_crai -> [meta, normal_cram, normal_crai, tumor_cram, tumor_crai, [], []] }

        BAM_VARIANT_CALLING_SOMATIC_STRELKA(
            cram_strelka,
            dict,
            fasta.map { _meta, fasta_ -> [fasta_] },
            fasta_fai.map { _meta, fasta_fai_ -> [fasta_fai_] },
            intervals_bed_gz_tbi,
        )

        vcf_strelka = BAM_VARIANT_CALLING_SOMATIC_STRELKA.out.vcf
        versions = versions.mix(BAM_VARIANT_CALLING_SOMATIC_STRELKA.out.versions)
    }

    // MSISENSOR2
    if (tools && tools.split(',').contains('msisensor2')) {
        // no need for models in tumor normal mode
        def models = []

        MSISENSOR2_MSI(bam.combine(intervals_bed_combined), msisensor2_scan, models)

        versions = versions.mix(MSISENSOR2_MSI.out.versions)
        out_msisensor2 = out_msisensor2.mix(MSISENSOR2_MSI.out.distribution)
        out_msisensor2 = out_msisensor2.mix(MSISENSOR2_MSI.out.somatic)
        out_msisensor2 = out_msisensor2.mix(MSISENSOR2_MSI.out.germline)
    }

    // MSISENSORPRO
    if (tools && tools.split(',').contains('msisensorpro')) {
        MSISENSORPRO_MSISOMATIC(cram.combine(intervals_bed_combined), fasta, msisensorpro_scan)

        versions = versions.mix(MSISENSORPRO_MSISOMATIC.out.versions)
        out_msisensorpro = out_msisensorpro.mix(MSISENSORPRO_MSISOMATIC.out.output_report)
    }

    // MuSE
    if (tools && tools.split(',').contains('muse')) {
        BAM_VARIANT_CALLING_SOMATIC_MUSE(
            bam.map { meta, normal_bam, normal_bai, _tumor_bam, _tumor_bai -> [meta, normal_bam, normal_bai] },
            bam.map { meta, _normal_bam, _normal_bai, tumor_bam, tumor_bai -> [meta, tumor_bam, tumor_bai] },
            fasta,
            dbsnp,
            dbsnp_tbi,
        )

        vcf_muse = BAM_VARIANT_CALLING_SOMATIC_MUSE.out.vcf
        versions = versions.mix(BAM_VARIANT_CALLING_SOMATIC_MUSE.out.versions)
    }

    // MUTECT2
    if (tools && tools.split(',').contains('mutect2')) {
        // joint_mutect2 mode needs different meta.map than regular mode
        //   we need to keep all fields and then remove on a per-tool-basis to ensure proper joining at the filtering step
        BAM_VARIANT_CALLING_SOMATIC_MUTECT2(
            cram.map { meta, normal_cram, normal_crai, tumor_cram, tumor_crai ->
                joint_mutect2
                    ? [meta + [id: meta.patient], [normal_cram, tumor_cram], [normal_crai, tumor_crai]]
                    : [meta, [normal_cram, tumor_cram], [normal_crai, tumor_crai]]
            },
            fasta,
            fasta_fai,
            dict,
            germline_resource,
            germline_resource_tbi,
            panel_of_normals,
            panel_of_normals_tbi,
            intervals,
            joint_mutect2,
        )

        vcf_mutect2 = BAM_VARIANT_CALLING_SOMATIC_MUTECT2.out.vcf_filtered
        versions = versions.mix(BAM_VARIANT_CALLING_SOMATIC_MUTECT2.out.versions)
    }

    // TNSCOPE
    if (tools && tools.split(',').contains('sentieon_tnscope')) {

        BAM_VARIANT_CALLING_SOMATIC_TNSCOPE(
            cram.map { meta, normal_cram, normal_crai, tumor_cram, tumor_crai ->
                [meta, [normal_cram, tumor_cram], [normal_crai, tumor_crai]]
            },
            fasta,
            fasta_fai,
            dict,
            germline_resource,
            germline_resource_tbi,
            panel_of_normals,
            panel_of_normals_tbi,
            intervals,
        )

        vcf_tnscope = BAM_VARIANT_CALLING_SOMATIC_TNSCOPE.out.vcf
        versions = versions.mix(BAM_VARIANT_CALLING_SOMATIC_TNSCOPE.out.versions)
    }

    // TIDDIT
    if (tools && tools.split(',').contains('tiddit')) {
        BAM_VARIANT_CALLING_SOMATIC_TIDDIT(
            cram.map { meta, normal_cram, normal_crai, _tumor_cram, _tumor_crai -> [meta, normal_cram, normal_crai] },
            cram.map { meta, _normal_cram, _normal_crai, tumor_cram, tumor_crai -> [meta, tumor_cram, tumor_crai] },
            fasta,
            bwa,
        )

        vcf_tiddit = BAM_VARIANT_CALLING_SOMATIC_TIDDIT.out.vcf
        versions = versions.mix(BAM_VARIANT_CALLING_SOMATIC_TIDDIT.out.versions)
    }

    vcf_all = Channel.empty()
        .mix(
            vcf_freebayes,
            vcf_manta,
            vcf_muse,
            vcf_mutect2,
            vcf_strelka,
            vcf_tiddit,
            vcf_tnscope,
        )

    emit:
    out_indexcov
    out_msisensor2
    out_msisensorpro
    vcf_all
    vcf_freebayes
    vcf_manta
    vcf_muse
    vcf_mutect2
    vcf_strelka
    vcf_tiddit
    vcf_tnscope
    versions
}
