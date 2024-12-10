//
// GERMLINE VARIANT CALLING
//

include { BAM_JOINT_CALLING_GERMLINE_GATK                                              } from '../bam_joint_calling_germline_gatk/main'
include { BAM_JOINT_CALLING_GERMLINE_SENTIEON                                          } from '../bam_joint_calling_germline_sentieon/main'
include { BAM_VARIANT_CALLING_CNVKIT                                                   } from '../bam_variant_calling_cnvkit/main'
include { BAM_VARIANT_CALLING_DEEPVARIANT                                              } from '../bam_variant_calling_deepvariant/main'
include { BAM_VARIANT_CALLING_FREEBAYES                                                } from '../bam_variant_calling_freebayes/main'
include { BAM_VARIANT_CALLING_GERMLINE_MANTA                                           } from '../bam_variant_calling_germline_manta/main'
include { BAM_VARIANT_CALLING_HAPLOTYPECALLER                                          } from '../bam_variant_calling_haplotypecaller/main'
include { BAM_VARIANT_CALLING_INDEXCOV                                                 } from '../bam_variant_calling_indexcov/main'
include { BAM_VARIANT_CALLING_SENTIEON_DNASCOPE                                        } from '../bam_variant_calling_sentieon_dnascope/main'
include { BAM_VARIANT_CALLING_SENTIEON_HAPLOTYPER                                      } from '../bam_variant_calling_sentieon_haplotyper/main'
include { BAM_VARIANT_CALLING_MPILEUP                                                  } from '../bam_variant_calling_mpileup/main'
include { BAM_VARIANT_CALLING_SINGLE_STRELKA                                           } from '../bam_variant_calling_single_strelka/main'
include { BAM_VARIANT_CALLING_SINGLE_TIDDIT                                            } from '../bam_variant_calling_single_tiddit/main'
include { SENTIEON_DNAMODELAPPLY                                                       } from '../../../modules/nf-core/sentieon/dnamodelapply/main'
include { VCF_VARIANT_FILTERING_GATK                                                   } from '../vcf_variant_filtering_gatk/main'
include { VCF_VARIANT_FILTERING_GATK as SENTIEON_HAPLOTYPER_VCF_VARIANT_FILTERING_GATK } from '../vcf_variant_filtering_gatk/main'

workflow BAM_VARIANT_CALLING_GERMLINE_ALL {
    take:
    tools                             // Mandatory, list of tools to apply
    skip_tools                        // Mandatory, list of tools to skip
    cram                              // channel: [mandatory] meta, cram
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
    skip_haplotypecaller_filter       // boolean: [mandatory] [default: false] whether to filter haplotypecaller single sample vcfs
    sentieon_haplotyper_emit_mode     // channel: [mandatory] value channel with string
    sentieon_dnascope_emit_mode       // channel: [mandatory] value channel with string
    sentieon_dnascope_pcr_indel_model // channel: [mandatory] value channel with string
    sentieon_dnascope_model           // channel: [mandatory] value channel with string

    main:
    versions = Channel.empty()

    //TODO: Temporary until the if's can be removed and printing to terminal is prevented with "when" in the modules.config
    gvcf_sentieon_dnascope   = Channel.empty()
    gvcf_sentieon_haplotyper = Channel.empty()

    out_indexcov             = Channel.empty()
    vcf_deepvariant          = Channel.empty()
    vcf_freebayes            = Channel.empty()
    vcf_haplotypecaller      = Channel.empty()
    vcf_manta                = Channel.empty()
    vcf_mpileup              = Channel.empty()
    vcf_sentieon_dnascope    = Channel.empty()
    vcf_sentieon_haplotyper  = Channel.empty()
    vcf_strelka              = Channel.empty()
    vcf_tiddit               = Channel.empty()

    // BCFTOOLS MPILEUP
    if (tools.split(',').contains('mpileup')) {
        BAM_VARIANT_CALLING_MPILEUP(
            cram,
            dict,
            fasta,
            intervals
        )
        vcf_mpileup = BAM_VARIANT_CALLING_MPILEUP.out.vcf
        versions = versions.mix(BAM_VARIANT_CALLING_MPILEUP.out.versions)
    }

    // CNVKIT
    if (tools.split(',').contains('cnvkit')) {
        BAM_VARIANT_CALLING_CNVKIT(
            // Remap channel to match module/subworkflow
            cram.map{ meta, cram, crai -> [ meta, [], cram ] },
            fasta,
            fasta_fai,
            intervals_bed_combined.map{ it -> [[id:it[0].baseName], it] },
            params.cnvkit_reference ? cnvkit_reference.map{ it -> [[id:it[0].baseName], it] } : [[:],[]]
        )
        versions = versions.mix(BAM_VARIANT_CALLING_CNVKIT.out.versions)
    }

    // DEEPVARIANT
    if (tools.split(',').contains('deepvariant')) {
        BAM_VARIANT_CALLING_DEEPVARIANT(
            cram,
            dict,
            fasta,
            fasta_fai,
            intervals
        )

        vcf_deepvariant = BAM_VARIANT_CALLING_DEEPVARIANT.out.vcf
        versions = versions.mix(BAM_VARIANT_CALLING_DEEPVARIANT.out.versions)
    }

    // FREEBAYES
    if (tools.split(',').contains('freebayes')) {
        // Input channel is remapped to match input of module/subworkflow
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

    // HAPLOTYPECALLER
    if (tools.split(',').contains('haplotypecaller')) {
        BAM_VARIANT_CALLING_HAPLOTYPECALLER(
            cram,
            fasta,
            fasta_fai,
            dict,
            dbsnp.map{it -> [[:], it]},
            dbsnp_tbi.map{it -> [[:], it]},
            intervals)

        vcf_haplotypecaller = BAM_VARIANT_CALLING_HAPLOTYPECALLER.out.vcf
        tbi_haplotypecaller = BAM_VARIANT_CALLING_HAPLOTYPECALLER.out.tbi

        versions = versions.mix(BAM_VARIANT_CALLING_HAPLOTYPECALLER.out.versions)

        if (joint_germline) {
            BAM_JOINT_CALLING_GERMLINE_GATK(
                BAM_VARIANT_CALLING_HAPLOTYPECALLER.out.gvcf_tbi_intervals,
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
                known_snps_vqsr)

            vcf_haplotypecaller = BAM_JOINT_CALLING_GERMLINE_GATK.out.genotype_vcf
            versions = versions.mix(BAM_JOINT_CALLING_GERMLINE_GATK.out.versions)
        } else {

            // If single sample track, check if filtering should be done
            if (!skip_haplotypecaller_filter) {

                VCF_VARIANT_FILTERING_GATK(
                    vcf_haplotypecaller.join(tbi_haplotypecaller, failOnDuplicate: true, failOnMismatch: true),
                    fasta.map{ meta, fasta -> [ fasta ] },
                    fasta_fai.map{ meta, fasta_fai -> [ fasta_fai ] },
                    dict.map{ meta, dict -> [ dict ] },
                    intervals_bed_combined_haplotypec,
                    known_sites_indels.concat(known_sites_snps).flatten().unique().collect(),
                    known_sites_indels_tbi.concat(known_sites_snps_tbi).flatten().unique().collect())

                vcf_haplotypecaller = VCF_VARIANT_FILTERING_GATK.out.filtered_vcf

                versions = versions.mix(VCF_VARIANT_FILTERING_GATK.out.versions)
            }
        }
    }

    // MANTA
    if (tools.split(',').contains('manta')) {
        BAM_VARIANT_CALLING_GERMLINE_MANTA (
            cram,
            fasta,
            fasta_fai,
            intervals_bed_gz_tbi_combined
        )

        vcf_manta = BAM_VARIANT_CALLING_GERMLINE_MANTA.out.vcf
        versions = versions.mix(BAM_VARIANT_CALLING_GERMLINE_MANTA.out.versions)
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

    // SENTIEON DNASCOPE
    if (tools.split(',').contains('sentieon_dnascope')) {
        BAM_VARIANT_CALLING_SENTIEON_DNASCOPE(
            cram,
            fasta,
            fasta_fai,
            dict,
            dbsnp,
            dbsnp_tbi,
            dbsnp_vqsr,
            intervals,
            joint_germline,
            sentieon_dnascope_emit_mode,
            sentieon_dnascope_pcr_indel_model,
            sentieon_dnascope_model)

        versions = versions.mix(BAM_VARIANT_CALLING_SENTIEON_DNASCOPE.out.versions)

        vcf_sentieon_dnascope      = BAM_VARIANT_CALLING_SENTIEON_DNASCOPE.out.vcf
        vcf_tbi_sentieon_dnascope  = BAM_VARIANT_CALLING_SENTIEON_DNASCOPE.out.vcf_tbi
        gvcf_sentieon_dnascope     = BAM_VARIANT_CALLING_SENTIEON_DNASCOPE.out.gvcf
        gvcf_tbi_sentieon_dnascope = BAM_VARIANT_CALLING_SENTIEON_DNASCOPE.out.gvcf_tbi

        if (joint_germline) {
            BAM_JOINT_CALLING_GERMLINE_SENTIEON(
                BAM_VARIANT_CALLING_SENTIEON_DNASCOPE.out.genotype_intervals,
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
                'sentieon_dnascope')

            vcf_sentieon_dnascope = BAM_JOINT_CALLING_GERMLINE_SENTIEON.out.genotype_vcf
            versions = versions.mix(BAM_JOINT_CALLING_GERMLINE_SENTIEON.out.versions)
        } else {
            // If single sample track, check if filtering should be done
            if (!(skip_tools && skip_tools.split(',').contains('dnascope_filter'))) {

                SENTIEON_DNAMODELAPPLY(
                    vcf_sentieon_dnascope.join(vcf_tbi_sentieon_dnascope, failOnDuplicate: true, failOnMismatch: true),
                    fasta,
                    fasta_fai,
                    sentieon_dnascope_model.map{ model -> [ [ id:model.baseName ], model ] })

                vcf_sentieon_dnascope = SENTIEON_DNAMODELAPPLY.out.vcf
                versions = versions.mix(SENTIEON_DNAMODELAPPLY.out.versions)

            }

        }
    }

    // SENTIEON HAPLOTYPER
    if (tools.split(',').contains('sentieon_haplotyper')) {
        BAM_VARIANT_CALLING_SENTIEON_HAPLOTYPER(
            cram,
            fasta,
            fasta_fai,
            dict,
            dbsnp,
            dbsnp_tbi,
            dbsnp_vqsr,
            intervals,
            joint_germline,
            sentieon_haplotyper_emit_mode)

        versions = versions.mix(BAM_VARIANT_CALLING_SENTIEON_HAPLOTYPER.out.versions)

        vcf_sentieon_haplotyper      = BAM_VARIANT_CALLING_SENTIEON_HAPLOTYPER.out.vcf
        vcf_tbi_sentieon_haplotyper  = BAM_VARIANT_CALLING_SENTIEON_HAPLOTYPER.out.vcf_tbi
        gvcf_sentieon_haplotyper     = BAM_VARIANT_CALLING_SENTIEON_HAPLOTYPER.out.gvcf
        gvcf_tbi_sentieon_haplotyper = BAM_VARIANT_CALLING_SENTIEON_HAPLOTYPER.out.gvcf_tbi

        if (joint_germline) {
            BAM_JOINT_CALLING_GERMLINE_SENTIEON(
                BAM_VARIANT_CALLING_SENTIEON_HAPLOTYPER.out.genotype_intervals,
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
                'sentieon_haplotyper')

            vcf_sentieon_haplotyper = BAM_JOINT_CALLING_GERMLINE_SENTIEON.out.genotype_vcf
            versions = versions.mix(BAM_JOINT_CALLING_GERMLINE_SENTIEON.out.versions)
        } else {

            // If single sample track, check if filtering should be done
            if (!(skip_tools && skip_tools.split(',').contains('haplotyper_filter'))) {

                SENTIEON_HAPLOTYPER_VCF_VARIANT_FILTERING_GATK(
                    vcf_sentieon_haplotyper.join(vcf_tbi_sentieon_haplotyper, failOnDuplicate: true, failOnMismatch: true),
                    fasta.map{ meta, it -> [ it ] },
                    fasta_fai.map{ meta, it -> [ it ] },
                    dict.map{ meta, dict -> [ dict ] },
                    intervals_bed_combined_haplotypec,
                    known_sites_indels.concat(known_sites_snps).flatten().unique().collect(),
                    known_sites_indels_tbi.concat(known_sites_snps_tbi).flatten().unique().collect())

                vcf_sentieon_haplotyper = SENTIEON_HAPLOTYPER_VCF_VARIANT_FILTERING_GATK.out.filtered_vcf

                versions = versions.mix(SENTIEON_HAPLOTYPER_VCF_VARIANT_FILTERING_GATK.out.versions)
            }
        }
    }

    // STRELKA

    if (tools.split(',').contains('strelka')) {

        BAM_VARIANT_CALLING_SINGLE_STRELKA(
            cram,
            dict,
            fasta.map{ meta, fasta -> [ fasta ] },
            fasta_fai.map{ meta, fasta_fai -> [ fasta_fai ] },
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
            fasta,
            bwa
        )

        vcf_tiddit = BAM_VARIANT_CALLING_SINGLE_TIDDIT.out.vcf
        versions = versions.mix(BAM_VARIANT_CALLING_SINGLE_TIDDIT.out.versions)
    }

    vcf_all = Channel.empty().mix(
        vcf_deepvariant,
        vcf_freebayes,
        vcf_sentieon_dnascope,
        vcf_haplotypecaller,
        vcf_manta,
        vcf_mpileup,
        vcf_sentieon_haplotyper,
        vcf_strelka,
        vcf_tiddit
    )

    emit:
    gvcf_sentieon_dnascope
    gvcf_sentieon_haplotyper
    out_indexcov
    vcf_all
    vcf_deepvariant
    vcf_freebayes
    vcf_haplotypecaller
    vcf_manta
    vcf_mpileup
    vcf_strelka
    vcf_sentieon_dnascope
    vcf_sentieon_haplotyper
    vcf_tiddit

    versions
}
