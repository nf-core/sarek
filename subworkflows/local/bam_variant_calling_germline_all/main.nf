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
    skip_tools                        // Mandatory, list of tools to skip
    cram                              // channel: [mandatory] cram
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
    versions = Channel.empty()

    //TODO: Temporary until the if's can be removed and printing to terminal is prevented with "when" in the modules.config
    vcf_deepvariant     = Channel.empty()
    vcf_freebayes       = Channel.empty()
    vcf_haplotypecaller = Channel.empty()
    vcf_manta           = Channel.empty()
    vcf_mpileup         = Channel.empty()
    vcf_strelka         = Channel.empty()
    vcf_tiddit          = Channel.empty()

    // Remap channel with intervals
    cram_intervals = cram
        .combine(intervals).map{ meta, cram, crai, intervals, num_intervals ->
            // If no interval file provided (0) then add empty list
            [ meta.subMap('data_type', 'patient', 'sample', 'sex', 'status')
                + [ id:meta.sample, num_intervals:num_intervals ],
                cram, crai, (num_intervals == 0 ? [] : intervals) ]
        }

    // Remap channel with gzipped intervals + indexes
    cram_intervals_gz_tbi = cram
        .combine(intervals_bed_gz_tbi).map{ meta, cram, crai, bed_tbi, num_intervals ->
            // If no interval file provided (0) then add empty list
            [ meta.subMap('data_type', 'patient', 'sample', 'sex', 'status')
                + [ id:meta.sample, num_intervals:num_intervals ],
                cram, crai, (num_intervals == 0 ? [] : bed_tbi[0]), (num_intervals == 0 ? [] : bed_tbi[1])]
        }

    if (tools.split(',').contains('mpileup')) {
        // Input channel is remapped to match input of module/subworkflow
        BAM_VARIANT_CALLING_MPILEUP(
            cram_intervals.map{ meta, cram, crai, intervals -> [ meta, cram, intervals ] },
            fasta,
            dict
        )

        vcf_mpileup = BAM_VARIANT_CALLING_MPILEUP.out.vcf
        versions = versions.mix(BAM_VARIANT_CALLING_MPILEUP.out.versions)
    }

    // CNVKIT

    if (tools.split(',').contains('cnvkit')) {
        // Input channel is remapped to match input of module/subworkflow
        BAM_VARIANT_CALLING_CNVKIT(
            cram.map{ meta, cram, crai -> [ meta, [], cram ] },
            fasta,
            fasta_fai,
            intervals_bed_combined,
            []
        )

        versions = versions.mix(BAM_VARIANT_CALLING_CNVKIT.out.versions)
    }

    // DEEPVARIANT
    if (tools.split(',').contains('deepvariant')) {
        BAM_VARIANT_CALLING_DEEPVARIANT(
            cram_intervals,
            dict,
            fasta,
            fasta_fai
        )

        vcf_deepvariant = BAM_VARIANT_CALLING_DEEPVARIANT.out.vcf
        versions = versions.mix(BAM_VARIANT_CALLING_DEEPVARIANT.out.versions)
    }

    // FREEBAYES
    if (tools.split(',').contains('freebayes')) {
        // Input channel is remapped to match input of module/subworkflow
        BAM_VARIANT_CALLING_FREEBAYES(
            cram_intervals.map{ meta, cram, crai, intervals -> [ meta, cram, crai, [], [], intervals ] },
            dict,
            fasta,
            fasta_fai
        )

        vcf_freebayes = BAM_VARIANT_CALLING_FREEBAYES.out.vcf
        versions = versions.mix(BAM_VARIANT_CALLING_FREEBAYES.out.versions)
    }

    // HAPLOTYPECALLER
    if (tools.split(',').contains('haplotypecaller')) {
        // Input channel is remapped to match input of module/subworkflow
        cram_intervals_haplotypecaller = cram_intervals
            .map{ meta, cram, crai, intervals ->
                [ (params.joint_germline ?
                    meta.subMap('data_type', 'num_intervals', 'patient', 'sample', 'sex', 'status')
                        + [ id:meta.sample, intervals_name:(meta.num_intervals == 0 ? "no_interval" : intervals.simpleName) ] :
                    meta),
                cram, crai, intervals, [] ]
        }

        BAM_VARIANT_CALLING_HAPLOTYPECALLER(
            cram_intervals_haplotypecaller,
            fasta,
            fasta_fai,
            dict,
            dbsnp,
            dbsnp_tbi,
            known_sites_indels,
            known_sites_indels_tbi,
            known_sites_snps,
            known_sites_snps_tbi,
            intervals_bed_combined_haplotypec,
            (skip_tools && skip_tools.split(',').contains('haplotypecaller_filter')))

        vcf_haplotypecaller = BAM_VARIANT_CALLING_HAPLOTYPECALLER.out.vcf
        versions = versions.mix(BAM_VARIANT_CALLING_HAPLOTYPECALLER.out.versions)
    }

    // MANTA
    if (tools.split(',').contains('manta')) {
        BAM_VARIANT_CALLING_GERMLINE_MANTA (
            cram_intervals_gz_tbi,
            dict,
            fasta,
            fasta_fai
        )

        vcf_manta = BAM_VARIANT_CALLING_GERMLINE_MANTA.out.vcf
        versions = versions.mix(BAM_VARIANT_CALLING_GERMLINE_MANTA.out.versions)
    }

    // STRELKA
    if (tools.split(',').contains('strelka')) {
        BAM_VARIANT_CALLING_SINGLE_STRELKA(
            cram_intervals_gz_tbi,
            dict,
            fasta,
            fasta_fai
        )

        vcf_strelka = BAM_VARIANT_CALLING_SINGLE_STRELKA.out.vcf
        versions = versions.mix(BAM_VARIANT_CALLING_SINGLE_STRELKA.out.versions)
    }

    //TIDDIT
    if (tools.split(',').contains('tiddit')) {
        BAM_VARIANT_CALLING_SINGLE_TIDDIT(
            cram,
            fasta.map{ it -> [[id:it[0].baseName], it] },
            bwa
        )

        vcf_tiddit = BAM_VARIANT_CALLING_SINGLE_TIDDIT.out.vcf
        versions = versions.mix(BAM_VARIANT_CALLING_SINGLE_TIDDIT.out.versions)
    }

    emit:
    vcf_deepvariant
    vcf_freebayes
    vcf_haplotypecaller
    vcf_manta
    vcf_mpileup
    vcf_strelka
    vcf_tiddit

    versions
}
