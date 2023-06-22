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
include { BAM_VARIANT_CALLING_TUMOR_ONLY_MUTECT2_MULTI_SAMPLE   } from '../bam_variant_calling_tumor_only_mutect2_ms/main'

workflow BAM_VARIANT_CALLING_TUMOR_ONLY_ALL {
    take:
    tools                         // Mandatory, list of tools to apply
    cram                          // channel: [mandatory] cram
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
    intervals                     // channel: [mandatory] [ intervals, num_intervals ] or [ [], 0 ] if no intervals
    intervals_bed_gz_tbi          // channel: [mandatory] [ interval.bed.gz, interval.bed.gz.tbi, num_intervals ] or [ [], [], 0 ] if no intervals
    intervals_bed_combined        // channel: [mandatory] intervals/target regions in one file unzipped
    intervals_bed_gz_tbi_combined // channel: [mandatory] intervals/target regions in one file zipped
    mappability
    panel_of_normals              // channel: [optional]  panel_of_normals
    panel_of_normals_tbi          // channel: [optional]  panel_of_normals_tbi

    main:
    versions = Channel.empty()

    //TODO: Temporary until the if's can be removed and printing to terminal is prevented with "when" in the modules.config
    vcf_freebayes   = Channel.empty()
    vcf_manta       = Channel.empty()
    vcf_mutect2     = Channel.empty()
    vcf_mutect2_ms  = Channel.empty()
    vcf_strelka     = Channel.empty()
    vcf_tiddit      = Channel.empty()

    // MPILEUP
    if (tools.split(',').contains('mpileup') || tools.split(',').contains('controlfreec')) {
        BAM_VARIANT_CALLING_MPILEUP(
            cram,
            dict,
            fasta,
            intervals
        )

        versions = versions.mix(BAM_VARIANT_CALLING_MPILEUP.out.versions)
    }

    // CONTROLFREEC (depends on MPILEUP)
    if (tools.split(',').contains('controlfreec')) {
        length_file = cf_chrom_len ?: fasta_fai

        BAM_VARIANT_CALLING_TUMOR_ONLY_CONTROLFREEC(
            // Remap channel to match module/subworkflow
            BAM_VARIANT_CALLING_MPILEUP.out.mpileup.map{ meta, pileup_tumor -> [ meta, [], pileup_tumor, [], [], [], [] ] },
            fasta,
            length_file,
            dbsnp,
            dbsnp_tbi,
            chr_files,
            mappability,
            intervals_bed_combined
        )

        versions = versions.mix(BAM_VARIANT_CALLING_TUMOR_ONLY_CONTROLFREEC.out.versions)
    }

    // CNVKIT
    if (tools.split(',').contains('cnvkit')) {
        BAM_VARIANT_CALLING_CNVKIT (
            // Remap channel to match module/subworkflow
            cram.map{ meta, cram, crai -> [ meta, cram, [] ] },
            fasta,
            fasta_fai,
            [],
            cnvkit_reference
        )

        versions = versions.mix(BAM_VARIANT_CALLING_CNVKIT.out.versions)
    }

    // FREEBAYES
    if (tools.split(',').contains('freebayes')) {
        BAM_VARIANT_CALLING_FREEBAYES(
            // Remap channel to match module/subworkflow
            cram.map{ meta, cram, crai -> [ meta, cram, crai, [], [] ] },
            dict,
            fasta,
            fasta_fai,
            intervals
        )

        vcf_freebayes = BAM_VARIANT_CALLING_FREEBAYES.out.vcf
        versions = versions.mix(BAM_VARIANT_CALLING_FREEBAYES.out.versions)
    }

    // MUTECT2
    if (tools.split(',').contains('mutect2')) {
        if (params.mutect2_multi_sample) {
            BAM_VARIANT_CALLING_TUMOR_ONLY_MUTECT2_MULTI_SAMPLE(
                cram,
                // Remap channel to match module/subworkflow
                fasta.map{ it -> [ [ id:'fasta' ], it ] },
                // Remap channel to match module/subworkflow
                fasta_fai.map{ it -> [ [ id:'fasta_fai' ], it ] },
                dict,
                germline_resource,
                germline_resource_tbi,
                panel_of_normals,
                panel_of_normals_tbi,
                intervals
            )
            vcf_mutect2_ms = BAM_VARIANT_CALLING_TUMOR_ONLY_MUTECT2_MULTI_SAMPLE.out.filtered_vcf
            versions = versions.mix(BAM_VARIANT_CALLING_TUMOR_ONLY_MUTECT2_MULTI_SAMPLE.out.versions)
        }
        else {
            BAM_VARIANT_CALLING_TUMOR_ONLY_MUTECT2(
                cram,
                // Remap channel to match module/subworkflow
                fasta.map{ it -> [ [ id:'fasta' ], it ] },
                // Remap channel to match module/subworkflow
                fasta_fai.map{ it -> [ [ id:'fasta_fai' ], it ] },
                dict,
                germline_resource,
                germline_resource_tbi,
                panel_of_normals,
                panel_of_normals_tbi,
                intervals
            )

            vcf_mutect2 = BAM_VARIANT_CALLING_TUMOR_ONLY_MUTECT2.out.vcf_filtered
            versions = versions.mix(BAM_VARIANT_CALLING_TUMOR_ONLY_MUTECT2.out.versions)
        }
    }

    // MANTA
    if (tools.split(',').contains('manta')) {
        BAM_VARIANT_CALLING_TUMOR_ONLY_MANTA(
            cram,
            // Remap channel to match module/subworkflow
            dict.map{ it -> [ [ id:'dict' ], it ] },
            fasta,
            fasta_fai,
            intervals_bed_gz_tbi_combined

        )

        vcf_manta = BAM_VARIANT_CALLING_TUMOR_ONLY_MANTA.out.vcf
        versions = versions.mix(BAM_VARIANT_CALLING_TUMOR_ONLY_MANTA.out.versions)
    }

    // STRELKA
    if (tools.split(',').contains('strelka')) {
        BAM_VARIANT_CALLING_SINGLE_STRELKA(
            cram,
            dict,
            fasta,
            fasta_fai,
            intervals_bed_gz_tbi
        )

        vcf_strelka = BAM_VARIANT_CALLING_SINGLE_STRELKA.out.vcf
        versions = versions.mix(BAM_VARIANT_CALLING_SINGLE_STRELKA.out.versions)
    }

    // TIDDIT
    if (tools.split(',').contains('tiddit')) {
        BAM_VARIANT_CALLING_SINGLE_TIDDIT(
            cram,
            // Remap channel to match module/subworkflow
            fasta.map{ it -> [ [ id:'fasta' ], it ] },
            bwa
        )

        vcf_tiddit = BAM_VARIANT_CALLING_SINGLE_TIDDIT.out.vcf
        versions = versions.mix(BAM_VARIANT_CALLING_SINGLE_TIDDIT.out.versions)
    }

    vcf_all = Channel.empty().mix(
        vcf_freebayes,
        vcf_manta,
        vcf_mutect2,
        vcf_mutect2_ms,
        vcf_strelka,
        vcf_tiddit
    )

    emit:
    vcf_all
    vcf_freebayes
    vcf_manta
    vcf_mutect2
    vcf_mutect2_ms
    vcf_strelka
    vcf_tiddit

    versions = versions
}
