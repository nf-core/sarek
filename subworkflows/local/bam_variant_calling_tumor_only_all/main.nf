//
// TUMOR ONLY VARIANT CALLING
// Should be only run on patients without normal sample
//

include { BAM_VARIANT_CALLING_CNVKIT                  } from '../bam_variant_calling_cnvkit/main'
include { BAM_VARIANT_CALLING_FREEBAYES               } from '../bam_variant_calling_freebayes/main'
include { BAM_VARIANT_CALLING_MPILEUP                 } from '../bam_variant_calling_mpileup/main'
include { BAM_VARIANT_CALLING_SINGLE_TIDDIT           } from '../bam_variant_calling_single_tiddit/main'
include { BAM_VARIANT_CALLING_TUMOR_ONLY_CONTROLFREEC } from '../bam_variant_calling_tumor_only_controlfreec/main'
include { BAM_VARIANT_CALLING_TUMOR_ONLY_MANTA        } from '../bam_variant_calling_tumor_only_manta/main'
include { BAM_VARIANT_CALLING_TUMOR_ONLY_MUTECT2      } from '../bam_variant_calling_tumor_only_mutect2/main'
include { BAM_VARIANT_CALLING_TUMOR_ONLY_LOFREQ       } from '../bam_variant_calling_tumor_only_lofreq/main'

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
    joint_mutect2                 // boolean: [mandatory] [default: false] run mutect2 in joint mode
    wes                           // boolean: [mandatory] [default: false] whether targeted data is processed

    main:
    versions = Channel.empty()

    //TODO: Temporary until the if's can be removed and printing to terminal is prevented with "when" in the modules.config
    vcf_freebayes   = Channel.empty()
    vcf_manta       = Channel.empty()
    vcf_mpileup     = Channel.empty()
    vcf_mutect2     = Channel.empty()
    vcf_tiddit      = Channel.empty()
    vcf_lofreq      = Channel.empty()

    // MPILEUP
    if (tools.split(',').contains('mpileup') || tools.split(',').contains('controlfreec')) {
        BAM_VARIANT_CALLING_MPILEUP(
            cram,
            dict,
            fasta,
            intervals
        )
        vcf_mpileup = BAM_VARIANT_CALLING_MPILEUP.out.vcf
        versions = versions.mix(BAM_VARIANT_CALLING_MPILEUP.out.versions)
    }

    // CONTROLFREEC (depends on MPILEUP)
    if (tools.split(',').contains('controlfreec')) {
        length_file            = cf_chrom_len ?: fasta_fai.map{ meta, fasta_fai -> [ fasta_fai ] }
        intervals_controlfreec = wes ? intervals_bed_combined : []

        BAM_VARIANT_CALLING_TUMOR_ONLY_CONTROLFREEC(
            // Remap channel to match module/subworkflow
            BAM_VARIANT_CALLING_MPILEUP.out.mpileup.map{ meta, pileup_tumor -> [ meta, [], pileup_tumor, [], [], [], [] ] },
            fasta.map{ meta, fasta -> [ fasta ] },
            length_file,
            dbsnp,
            dbsnp_tbi,
            chr_files,
            mappability,
            intervals_controlfreec
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
            [[id:"null"], []],
            cnvkit_reference.map{ it -> [[id:it[0].baseName], it] }
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
        BAM_VARIANT_CALLING_TUMOR_ONLY_MUTECT2(
            // Adjust meta.map to simplify joining channels
            cram.map{ meta, cram, crai ->
                joint_mutect2 ?
                //we need to keep all fields and then remove on a per-tool-basis to ensure proper joining at the filtering step
                [ meta - meta.subMap('data_type', 'status') + [ id:meta.patient ], cram, crai ] :
                [ meta - meta.subMap('data_type', 'status'), cram, crai ]
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

        vcf_mutect2 = BAM_VARIANT_CALLING_TUMOR_ONLY_MUTECT2.out.vcf_filtered
        versions = versions.mix(BAM_VARIANT_CALLING_TUMOR_ONLY_MUTECT2.out.versions)
    }

    //LOFREQ
    if (tools.split(',').contains('lofreq')) {
        BAM_VARIANT_CALLING_TUMOR_ONLY_LOFREQ(
            cram,
            fasta,
            fasta_fai,
            intervals,
            dict
        )
        vcf_lofreq = BAM_VARIANT_CALLING_TUMOR_ONLY_LOFREQ.out.vcf
        versions = versions.mix(BAM_VARIANT_CALLING_TUMOR_ONLY_LOFREQ.out.versions)
    }

    // MANTA
    if (tools.split(',').contains('manta')) {
        BAM_VARIANT_CALLING_TUMOR_ONLY_MANTA(
            cram,
            // Remap channel to match module/subworkflow
            fasta,
            fasta_fai,
            intervals_bed_gz_tbi_combined

        )

        vcf_manta = BAM_VARIANT_CALLING_TUMOR_ONLY_MANTA.out.vcf
        versions = versions.mix(BAM_VARIANT_CALLING_TUMOR_ONLY_MANTA.out.versions)
    }

    // TIDDIT
    if (tools.split(',').contains('tiddit')) {
        BAM_VARIANT_CALLING_SINGLE_TIDDIT(
            cram,
            fasta,
            bwa)

        vcf_tiddit = BAM_VARIANT_CALLING_SINGLE_TIDDIT.out.vcf
        versions = versions.mix(BAM_VARIANT_CALLING_SINGLE_TIDDIT.out.versions)
    }

    vcf_all = Channel.empty().mix(
        vcf_freebayes,
        vcf_lofreq,
        vcf_manta,
        vcf_mutect2,
        vcf_mpileup,
        vcf_tiddit
    )

    emit:
    vcf_all
    vcf_freebayes
    vcf_lofreq
    vcf_manta
    vcf_mpileup
    vcf_mutect2
    vcf_tiddit

    versions = versions
}
