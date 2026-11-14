//
// GERMLINE VARIANT CALLING
//

include { BAM_VARIANT_CALLING_CNVKIT                                                                                  } from '../bam_variant_calling_cnvkit/main'
include { BAM_VARIANT_CALLING_DEEPVARIANT                                                                             } from '../bam_variant_calling_deepvariant/main'
include { BAM_VARIANT_CALLING_PARABRICKS_DEEPVARIANT                                                                  } from '../bam_variant_calling_parabricks_deepvariant/main'
include { BAM_VARIANT_CALLING_FREEBAYES                                                                               } from '../bam_variant_calling_freebayes/main'
include { BAM_VARIANT_CALLING_GERMLINE_MANTA                                                                          } from '../bam_variant_calling_germline_manta/main'
include { BAM_VARIANT_CALLING_HAPLOTYPECALLER                                                                         } from '../bam_variant_calling_haplotypecaller/main'
include { BAM_VARIANT_CALLING_PARABRICKS_HAPLOTYPECALLER                                                              } from '../bam_variant_calling_parabricks_haplotypecaller/main'
include { BAM_VARIANT_CALLING_INDEXCOV                                                                                } from '../bam_variant_calling_indexcov/main'
include { BAM_VARIANT_CALLING_SENTIEON_DNASCOPE                                                                       } from '../bam_variant_calling_sentieon_dnascope/main'
include { BAM_VARIANT_CALLING_SENTIEON_HAPLOTYPER                                                                     } from '../bam_variant_calling_sentieon_haplotyper/main'
include { BAM_VARIANT_CALLING_MPILEUP                                                                                 } from '../bam_variant_calling_mpileup/main'
include { BAM_VARIANT_CALLING_SINGLE_STRELKA                                                                          } from '../bam_variant_calling_single_strelka/main'
include { BAM_VARIANT_CALLING_SINGLE_TIDDIT                                                                           } from '../bam_variant_calling_single_tiddit/main'
include { VCF_JOINT_CALLING_OR_FILTERING_HAPLOTYPECALLER                                                              } from '../vcf_joint_calling_or_filtering_haplotypecaller/main'
include { VCF_JOINT_CALLING_OR_FILTERING_HAPLOTYPECALLER as PARABRICKS_HAPLOTYPECALLER_VCF_JOINT_CALLING_OR_FILTERING } from '../vcf_joint_calling_or_filtering_haplotypecaller/main'
include { VCF_JOINT_CALLING_OR_FILTERING_SENTIEON                                                                     } from '../vcf_joint_calling_or_filtering_sentieon/main'
include { VCF_JOINT_CALLING_OR_FILTERING_SENTIEON as SENTIEON_HAPLOTYPER_VCF_JOINT_CALLING_OR_FILTERING               } from '../vcf_joint_calling_or_filtering_sentieon/main'

workflow BAM_VARIANT_CALLING_GERMLINE_ALL {
    take:
    tools                             // Mandatory, list of tools to apply
    skip_tools                        // Mandatory, list of tools to skip
    bam                               // channel: [mandatory] meta, bam
    cram                              // channel: [mandatory] meta, cram, crai
    bwa                               // channel: [mandatory] meta, bwa
    cnvkit_reference                  // channel: [optional] cnvkit reference
    dbsnp                             // channel: [mandatory] meta, dbsnp
    dbsnp_tbi                         // channel: [mandatory] dbsnp_tbi
    dbsnp_vqsr
    dict                              // channel: [mandatory] meta, dict
    fasta                             // channel: [mandatory] meta, fasta
    fasta_fai                         // channel: [mandatory] meta, fasta_fai
    intervals                         // channel: [mandatory] [ intervals, num_intervals ] or [ [], 0 ] if no intervals
    intervals_bed_combined            // channel: [mandatory] intervals/target regions in one file unzipped
    intervals_bed_gz_tbi_combined     // channel: [mandatory] intervals/target regions in one file zipped
    intervals_bed_combined_haplotypec // channel: [mandatory] intervals/target regions in one file unzipped, no_intervals.bed if no_intervals
    intervals_bed_gz_tbi              // channel: [mandatory] [ interval.bed.gz, interval.bed.gz.tbi, num_intervals ] or [ [], [], 0 ] if no intervals
    known_indels_vqsr
    known_sites_indels
    known_sites_indels_tbi
    known_sites_snps
    known_sites_snps_tbi
    known_snps_vqsr
    joint_germline                    // boolean: [mandatory] [default: false] joint calling of germline variants
    wes                               // boolean: [mandatory] [default: false] whether targeted data is processed
    skip_haplotypecaller_filter       // boolean: [mandatory] [default: false] whether to filter haplotypecaller single sample vcfs
    filter_resource_provided          // boolean: [mandatory] dbsnp, known indels or known snps given, needed by FilterVariantTranches
    sentieon_haplotyper_emit_mode     // channel: [mandatory] value channel with string
    sentieon_dnascope_emit_mode       // channel: [mandatory] value channel with string
    sentieon_dnascope_pcr_indel_model // channel: [mandatory] value channel with string
    sentieon_dnascope_model           // channel: [mandatory] value channel with string

    main:
    //TODO: Temporary until the if's can be removed and printing to terminal is prevented with "when" in the modules.config
    gvcf_parabricks_haplotypecaller     = channel.empty()
    gvcf_sentieon_dnascope              = channel.empty()
    gvcf_sentieon_haplotyper            = channel.empty()
    gvcf_tbi_parabricks_haplotypecaller = channel.empty()
    gvcf_tbi_sentieon_dnascope          = channel.empty()
    gvcf_tbi_sentieon_haplotyper        = channel.empty()

    out_indexcov                       = channel.empty()
    vcf_deepvariant                    = channel.empty()
    vcf_parabricks_deepvariant         = channel.empty()
    vcf_freebayes                      = channel.empty()
    vcf_haplotypecaller                = channel.empty()
    vcf_manta                          = channel.empty()
    vcf_mpileup                        = channel.empty()
    vcf_parabricks_haplotypecaller     = channel.empty()
    vcf_sentieon_dnascope              = channel.empty()
    vcf_sentieon_haplotyper            = channel.empty()
    vcf_strelka                        = channel.empty()
    vcf_tiddit                         = channel.empty()
    tbi_deepvariant                    = channel.empty()
    tbi_parabricks_deepvariant         = channel.empty()
    tbi_freebayes                      = channel.empty()
    tbi_haplotypecaller                = channel.empty()
    tbi_manta                          = channel.empty()
    tbi_mpileup                        = channel.empty()
    tbi_parabricks_haplotypecaller     = channel.empty()
    tbi_sentieon_dnascope              = channel.empty()
    tbi_sentieon_haplotyper            = channel.empty()
    tbi_strelka                        = channel.empty()
    tbi_tiddit                         = channel.empty()

    // BCFTOOLS MPILEUP
    if ('mpileup' in tools) {
        BAM_VARIANT_CALLING_MPILEUP(
            cram,
            dict,
            fasta,
            fasta_fai,
            intervals,
            tools,
        )
        vcf_mpileup = BAM_VARIANT_CALLING_MPILEUP.out.vcf
        tbi_mpileup = BAM_VARIANT_CALLING_MPILEUP.out.tbi
    }

    // CNVKIT
    if ('cnvkit' in tools) {
        BAM_VARIANT_CALLING_CNVKIT(
            // Remap channel to match module/subworkflow
            // Use the already-converted BAM (shared with somatic/tumor-only) instead of CRAM
            bam.map{ meta, bam_, _bai -> [ meta, [], bam_ ] },
            fasta,
            fasta_fai,
            intervals_bed_combined.map{_intervals -> _intervals ? [[id:_intervals[0].baseName], _intervals]: [[id:'no_intervals'], []]},
            params.cnvkit_reference ? cnvkit_reference.map{ reference -> [[id:reference[0].baseName], reference] } : [[:],[]]
        )
    }

    // DEEPVARIANT
    if ('deepvariant' in tools) {
        BAM_VARIANT_CALLING_DEEPVARIANT(
            cram,
            dict,
            fasta,
            fasta_fai,
            intervals
        )

        vcf_deepvariant = BAM_VARIANT_CALLING_DEEPVARIANT.out.vcf
        tbi_deepvariant = BAM_VARIANT_CALLING_DEEPVARIANT.out.tbi
    }

    // PARABRICKS DEEPVARIANT
    if ('parabricks_deepvariant' in tools) {
        BAM_VARIANT_CALLING_PARABRICKS_DEEPVARIANT(
            cram,
            fasta,
            intervals_bed_combined
        )

        vcf_parabricks_deepvariant = BAM_VARIANT_CALLING_PARABRICKS_DEEPVARIANT.out.vcf
        tbi_parabricks_deepvariant = BAM_VARIANT_CALLING_PARABRICKS_DEEPVARIANT.out.tbi
    }

    // FREEBAYES
    if ('freebayes' in tools) {
        // Input channel is remapped to match input of module/subworkflow
        BAM_VARIANT_CALLING_FREEBAYES(
            // Remap channel to match module/subworkflow
            cram.map{ meta, cram_, crai -> [ meta, cram_, crai, [], [] ] },
            dict,
            fasta,
            fasta_fai,
            intervals
        )

        vcf_freebayes = BAM_VARIANT_CALLING_FREEBAYES.out.vcf
        tbi_freebayes = BAM_VARIANT_CALLING_FREEBAYES.out.tbi
    }

    // HAPLOTYPECALLER
    if ('haplotypecaller' in tools) {
        BAM_VARIANT_CALLING_HAPLOTYPECALLER(
            cram,
            fasta,
            fasta_fai,
            dict,
            dbsnp.map{dbsnp_ -> [[:], dbsnp_]},
            dbsnp_tbi.map{dbsnp_tbi_ -> [[:], dbsnp_tbi_]},
            intervals)

        VCF_JOINT_CALLING_OR_FILTERING_HAPLOTYPECALLER(
            BAM_VARIANT_CALLING_HAPLOTYPECALLER.out.gvcf_tbi_intervals,
            BAM_VARIANT_CALLING_HAPLOTYPECALLER.out.vcf,
            BAM_VARIANT_CALLING_HAPLOTYPECALLER.out.tbi,
            fasta,
            fasta_fai,
            dict,
            dbsnp,
            dbsnp_tbi,
            dbsnp_vqsr,
            known_sites_indels,
            known_sites_indels_tbi,
            known_indels_vqsr,
            known_sites_snps,
            known_sites_snps_tbi,
            known_snps_vqsr,
            intervals_bed_combined_haplotypec,
            'haplotypecaller',
            joint_germline,
            skip_haplotypecaller_filter || !filter_resource_provided)

        vcf_haplotypecaller = VCF_JOINT_CALLING_OR_FILTERING_HAPLOTYPECALLER.out.vcf
        tbi_haplotypecaller = VCF_JOINT_CALLING_OR_FILTERING_HAPLOTYPECALLER.out.tbi
    }

    // PARABRICKS HAPLOTYPECALLER
    if ('parabricks_haplotypecaller' in tools) {
        BAM_VARIANT_CALLING_PARABRICKS_HAPLOTYPECALLER(
            cram,
            fasta,
            dbsnp,
            dbsnp_tbi,
            intervals,
            intervals_bed_combined,
            joint_germline
        )

        gvcf_parabricks_haplotypecaller     = BAM_VARIANT_CALLING_PARABRICKS_HAPLOTYPECALLER.out.gvcf
        gvcf_tbi_parabricks_haplotypecaller = BAM_VARIANT_CALLING_PARABRICKS_HAPLOTYPECALLER.out.gvcf_tbi

        PARABRICKS_HAPLOTYPECALLER_VCF_JOINT_CALLING_OR_FILTERING(
            BAM_VARIANT_CALLING_PARABRICKS_HAPLOTYPECALLER.out.gvcf_tbi_intervals,
            BAM_VARIANT_CALLING_PARABRICKS_HAPLOTYPECALLER.out.vcf,
            BAM_VARIANT_CALLING_PARABRICKS_HAPLOTYPECALLER.out.tbi,
            fasta,
            fasta_fai,
            dict,
            dbsnp,
            dbsnp_tbi,
            dbsnp_vqsr,
            known_sites_indels,
            known_sites_indels_tbi,
            known_indels_vqsr,
            known_sites_snps,
            known_sites_snps_tbi,
            known_snps_vqsr,
            intervals_bed_combined_haplotypec,
            'parabricks_haplotypecaller',
            joint_germline,
            'parabricks_haplotypecaller_filter' in skip_tools || !filter_resource_provided)

        vcf_parabricks_haplotypecaller = PARABRICKS_HAPLOTYPECALLER_VCF_JOINT_CALLING_OR_FILTERING.out.vcf
        tbi_parabricks_haplotypecaller = PARABRICKS_HAPLOTYPECALLER_VCF_JOINT_CALLING_OR_FILTERING.out.tbi
    }

    // MANTA
    if ('manta' in tools) {
        BAM_VARIANT_CALLING_GERMLINE_MANTA (
            cram,
            fasta,
            fasta_fai,
            intervals_bed_gz_tbi_combined
        )

        vcf_manta = BAM_VARIANT_CALLING_GERMLINE_MANTA.out.diploid_sv_vcf
        tbi_manta = BAM_VARIANT_CALLING_GERMLINE_MANTA.out.diploid_sv_vcf_tbi
    }

    // INDEXCOV, for WGS only
    if (wes == false && 'indexcov' in tools) {
        BAM_VARIANT_CALLING_INDEXCOV (
            cram,
            fasta,
            fasta_fai
        )

        out_indexcov = BAM_VARIANT_CALLING_INDEXCOV.out.out_indexcov
    }

    // SENTIEON DNASCOPE
    if ('sentieon_dnascope' in tools) {
        BAM_VARIANT_CALLING_SENTIEON_DNASCOPE(
            cram,
            fasta,
            fasta_fai,
            dict,
            dbsnp,
            dbsnp_tbi,
            intervals,
            joint_germline,
            sentieon_dnascope_emit_mode,
            sentieon_dnascope_pcr_indel_model,
            sentieon_dnascope_model)

        gvcf_sentieon_dnascope     = BAM_VARIANT_CALLING_SENTIEON_DNASCOPE.out.gvcf
        gvcf_tbi_sentieon_dnascope = BAM_VARIANT_CALLING_SENTIEON_DNASCOPE.out.gvcf_tbi

        VCF_JOINT_CALLING_OR_FILTERING_SENTIEON(
            BAM_VARIANT_CALLING_SENTIEON_DNASCOPE.out.genotype_intervals,
            BAM_VARIANT_CALLING_SENTIEON_DNASCOPE.out.vcf,
            BAM_VARIANT_CALLING_SENTIEON_DNASCOPE.out.vcf_tbi,
            fasta,
            fasta_fai,
            dict,
            dbsnp,
            dbsnp_tbi,
            dbsnp_vqsr,
            known_sites_indels,
            known_sites_indels_tbi,
            known_indels_vqsr,
            known_sites_snps,
            known_sites_snps_tbi,
            known_snps_vqsr,
            intervals_bed_combined_haplotypec,
            sentieon_dnascope_model,
            'sentieon_dnascope',
            joint_germline,
            'dnascope_filter' in skip_tools)

        vcf_sentieon_dnascope = VCF_JOINT_CALLING_OR_FILTERING_SENTIEON.out.vcf
        tbi_sentieon_dnascope = VCF_JOINT_CALLING_OR_FILTERING_SENTIEON.out.tbi
    }

    // SENTIEON HAPLOTYPER
    if ('sentieon_haplotyper' in tools) {
        BAM_VARIANT_CALLING_SENTIEON_HAPLOTYPER(
            cram,
            fasta,
            fasta_fai,
            dict,
            dbsnp,
            dbsnp_tbi,
            intervals,
            joint_germline,
            sentieon_haplotyper_emit_mode)

        gvcf_sentieon_haplotyper     = BAM_VARIANT_CALLING_SENTIEON_HAPLOTYPER.out.gvcf
        gvcf_tbi_sentieon_haplotyper = BAM_VARIANT_CALLING_SENTIEON_HAPLOTYPER.out.gvcf_tbi

        SENTIEON_HAPLOTYPER_VCF_JOINT_CALLING_OR_FILTERING(
            BAM_VARIANT_CALLING_SENTIEON_HAPLOTYPER.out.genotype_intervals,
            BAM_VARIANT_CALLING_SENTIEON_HAPLOTYPER.out.vcf,
            BAM_VARIANT_CALLING_SENTIEON_HAPLOTYPER.out.vcf_tbi,
            fasta,
            fasta_fai,
            dict,
            dbsnp,
            dbsnp_tbi,
            dbsnp_vqsr,
            known_sites_indels,
            known_sites_indels_tbi,
            known_indels_vqsr,
            known_sites_snps,
            known_sites_snps_tbi,
            known_snps_vqsr,
            intervals_bed_combined_haplotypec,
            sentieon_dnascope_model,
            'sentieon_haplotyper',
            joint_germline,
            'haplotyper_filter' in skip_tools || !filter_resource_provided)

        vcf_sentieon_haplotyper = SENTIEON_HAPLOTYPER_VCF_JOINT_CALLING_OR_FILTERING.out.vcf
        tbi_sentieon_haplotyper = SENTIEON_HAPLOTYPER_VCF_JOINT_CALLING_OR_FILTERING.out.tbi
    }

    // STRELKA

    if ('strelka' in tools) {

        BAM_VARIANT_CALLING_SINGLE_STRELKA(
            cram,
            dict,
            fasta.map{ _meta, fasta_ -> [ fasta_ ] },
            fasta_fai.map{ _meta, fasta_fai_ -> [ fasta_fai_ ] },
            intervals_bed_gz_tbi
        )

        vcf_strelka = BAM_VARIANT_CALLING_SINGLE_STRELKA.out.vcf
        tbi_strelka = BAM_VARIANT_CALLING_SINGLE_STRELKA.out.tbi
    }

    // TIDDIT
    if ('tiddit' in tools) {
        BAM_VARIANT_CALLING_SINGLE_TIDDIT(
            cram,
            // Remap channel to match module/subworkflow
            fasta,
            fasta_fai,
            bwa
        )

        vcf_tiddit = BAM_VARIANT_CALLING_SINGLE_TIDDIT.out.vcf
        tbi_tiddit = BAM_VARIANT_CALLING_SINGLE_TIDDIT.out.tbi
    }

    vcf_all = channel.empty().mix(
        vcf_deepvariant,
        vcf_parabricks_deepvariant,
        vcf_freebayes,
        vcf_haplotypecaller,
        vcf_manta,
        vcf_mpileup,
        vcf_parabricks_haplotypecaller,
        vcf_sentieon_dnascope,
        vcf_sentieon_haplotyper,
        vcf_strelka,
        vcf_tiddit
    )

    tbi_all = channel.empty().mix(
        tbi_deepvariant,
        tbi_parabricks_deepvariant,
        tbi_freebayes,
        tbi_haplotypecaller,
        tbi_manta,
        tbi_mpileup,
        tbi_parabricks_haplotypecaller,
        tbi_sentieon_dnascope,
        tbi_sentieon_haplotyper,
        tbi_strelka,
        tbi_tiddit
    )

    emit:
    gvcf_parabricks_haplotypecaller
    gvcf_sentieon_dnascope
    gvcf_sentieon_haplotyper
    gvcf_tbi_parabricks_haplotypecaller
    gvcf_tbi_sentieon_dnascope
    gvcf_tbi_sentieon_haplotyper
    out_indexcov
    vcf_all
    vcf_deepvariant
    vcf_parabricks_deepvariant
    vcf_freebayes
    vcf_haplotypecaller
    vcf_manta
    vcf_mpileup
    vcf_parabricks_haplotypecaller
    vcf_strelka
    vcf_sentieon_dnascope
    vcf_sentieon_haplotyper
    vcf_tiddit
    tbi_all
    tbi_deepvariant
    tbi_parabricks_deepvariant
    tbi_freebayes
    tbi_haplotypecaller
    tbi_manta
    tbi_mpileup
    tbi_parabricks_haplotypecaller
    tbi_sentieon_dnascope
    tbi_sentieon_haplotyper
    tbi_strelka
    tbi_tiddit
}
