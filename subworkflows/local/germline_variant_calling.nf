//
// GERMLINE VARIANT CALLING
//

include { RUN_CNVKIT_GERMLINE } from '../nf-core/variantcalling/cnvkit/germline/main.nf'
include { RUN_DEEPVARIANT     } from '../nf-core/variantcalling/deepvariant/main.nf'
include { RUN_FREEBAYES       } from '../nf-core/variantcalling/freebayes/main.nf'
include { RUN_HAPLOTYPECALLER } from '../nf-core/variantcalling/haplotypecaller/main.nf'
include { RUN_MANTA_GERMLINE  } from '../nf-core/variantcalling/manta/germline/main.nf'
include { RUN_MPILEUP         } from '../nf-core/variantcalling/mpileup/main'
include { RUN_STRELKA_SINGLE  } from '../nf-core/variantcalling/strelka/single/main.nf'
include { RUN_TIDDIT          } from '../nf-core/variantcalling/tiddit/main.nf'

workflow GERMLINE_VARIANT_CALLING {
    take:
        tools                         // Mandatory, list of tools to apply
        cram_recalibrated             // channel: [mandatory] cram
        bwa                           // channel: [mandatory] bwa
        dbsnp                         // channel: [mandatory] dbsnp
        dbsnp_tbi                     // channel: [mandatory] dbsnp_tbi
        dict                          // channel: [mandatory] dict
        fasta                         // channel: [mandatory] fasta
        fasta_fai                     // channel: [mandatory] fasta_fai
        intervals                     // channel: [mandatory] intervals/target regions
        intervals_bed_gz_tbi          // channel: [mandatory] intervals/target regions index zipped and indexed
        intervals_bed_combined        // channel: [mandatory] intervals/target regions in one file unzipped
        intervals_bed_combined_haplotypec // channel: [mandatory] intervals/target regions in one file unzipped, no_intervals.bed if no_intervals
        known_sites_indels
        known_sites_indels_tbi
        known_sites_snps
        known_sites_snps_tbi
        // joint_germline                // val: true/false on whether to run joint_germline calling, only works in combination with haplotypecaller at the moment

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
            intervals_new  = num_intervals == 0 ? []            : intervals
            intervals_name = num_intervals == 0 ? "no_interval" : intervals.simpleName

            [[patient:meta.patient, sample:meta.sample, sex:meta.sex, status:meta.status, id:meta.sample, data_type:meta.data_type, num_intervals:num_intervals, intervals_name:intervals_name],
            cram, crai, intervals_new]
        }

    // Remap channel with gzipped intervals + indexes
    cram_recalibrated_intervals_gz_tbi = cram_recalibrated.combine(intervals_bed_gz_tbi)
        .map{ meta, cram, crai, bed_tbi, num_intervals ->

            //If no interval file provided (0) then add empty list
            bed_new = num_intervals == 0 ? [] : bed_tbi[0]
            tbi_new = num_intervals == 0 ? [] : bed_tbi[1]

            [[patient:meta.patient, sample:meta.sample, sex:meta.sex, status:meta.status, id:meta.sample, data_type:meta.data_type, num_intervals:num_intervals],
            cram, crai, bed_new, tbi_new]
        }

    if(params.tools.contains('mpileup')){
        cram_intervals_no_index = cram_recalibrated_intervals
            .map { meta, cram, crai, intervals ->
                [meta, cram, intervals]
            }

        RUN_MPILEUP(cram_intervals_no_index,
                        fasta)
        mpileup_germline = RUN_MPILEUP.out.mpileup
        ch_versions = ch_versions.mix(RUN_MPILEUP.out.versions)
    }

    // CNVKIT

    if(tools.contains('cnvkit')){
        cram_recalibrated_cnvkit_germline = cram_recalibrated
            .map{ meta, cram, crai ->
                [meta, [], cram]
            }

        RUN_CNVKIT_GERMLINE(cram_recalibrated_cnvkit_germline,
                            fasta,
                            fasta_fai,
                            intervals_bed_combined,
                            [])
        ch_versions     = ch_versions.mix(RUN_CNVKIT_GERMLINE.out.versions)
    }

    // DEEPVARIANT
    if(tools.contains('deepvariant')){
        RUN_DEEPVARIANT(cram_recalibrated_intervals, dict, fasta, fasta_fai)

        deepvariant_vcf = Channel.empty().mix(RUN_DEEPVARIANT.out.deepvariant_vcf,RUN_DEEPVARIANT.out.deepvariant_gvcf)
        ch_versions     = ch_versions.mix(RUN_DEEPVARIANT.out.versions)
    }

    // FREEBAYES
    if (tools.contains('freebayes')){
        // Remap channel for Freebayes
        cram_recalibrated_intervals_freebayes = cram_recalibrated_intervals
            .map{ meta, cram, crai, intervals ->
                [meta, cram, crai, [], [], intervals]
            }
        RUN_FREEBAYES(cram_recalibrated_intervals_freebayes, dict, fasta, fasta_fai)

        freebayes_vcf   = RUN_FREEBAYES.out.freebayes_vcf
        ch_versions     = ch_versions.mix(RUN_FREEBAYES.out.versions)
    }

    // HAPLOTYPECALLER
    if (tools.contains('haplotypecaller')){
        RUN_HAPLOTYPECALLER(cram_recalibrated_intervals,
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

        haplotypecaller_vcf  = RUN_HAPLOTYPECALLER.out.filtered_vcf
        ch_versions          = ch_versions.mix(RUN_HAPLOTYPECALLER.out.versions)
    }

    // MANTA
    if (tools.contains('manta')){
        RUN_MANTA_GERMLINE (cram_recalibrated_intervals_gz_tbi,
                        dict,
                        fasta,
                        fasta_fai)

        manta_vcf   = RUN_MANTA_GERMLINE.out.manta_vcf
        ch_versions = ch_versions.mix(RUN_MANTA_GERMLINE.out.versions)
    }

    // STRELKA
    if (tools.contains('strelka')){
        RUN_STRELKA_SINGLE(cram_recalibrated_intervals_gz_tbi,
                dict,
                fasta,
                fasta_fai)

        strelka_vcf = RUN_STRELKA_SINGLE.out.strelka_vcf
        ch_versions = ch_versions.mix(RUN_STRELKA_SINGLE.out.versions)
    }

    //TIDDIT
    if (tools.contains('tiddit')){
        RUN_TIDDIT(cram_recalibrated,
                fasta,
                bwa)

        tiddit_vcf = RUN_TIDDIT.out.tiddit_vcf
        ch_versions = ch_versions.mix(RUN_TIDDIT.out.versions)
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
