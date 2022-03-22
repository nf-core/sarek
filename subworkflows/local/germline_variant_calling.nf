//
// GERMLINE VARIANT CALLING
//

include { DEEPVARIANT     } from './variantcalling/deepvariant.nf'
include { FREEBAYES       } from './variantcalling/freebayes.nf'
include { HAPLOTYPECALLER } from './variantcalling/haplotypecaller.nf'
include { MANTA_GERMLINE  } from './variantcalling/manta_germline.nf'
include { STRELKA_SINGLE  } from './variantcalling/strelka_single.nf'
//include { TIDDIT          } from './variantcalling/tiddit.nf'

workflow GERMLINE_VARIANT_CALLING {
    take:
        cram_recalibrated            // channel: [mandatory] cram
        dbsnp                        // channel: [mandatory] dbsnp
        dbsnp_tbi                    // channel: [mandatory] dbsnp_tbi
        dict                         // channel: [mandatory] dict
        fasta                        // channel: [mandatory] fasta
        fasta_fai                    // channel: [mandatory] fasta_fai
        intervals                    // channel: [mandatory] intervals/target regions
        intervals_bed_gz_tbi         // channel: [mandatory] intervals/target regions index zipped and indexed
        intervals_bed_combine_gz_tbi // channel: [mandatory] intervals/target regions index zipped and indexed in one file
        intervals_bed_combine_gz     // channel: [mandatory] intervals/target regions index zipped in one file
        num_intervals                // val: number of intervals that are used to parallelize exection, either based on capture kit or GATK recommended for WGS
        // joint_germline               // val: true/false on whether to run joint_germline calling, only works in combination with haplotypecaller at the moment

    main:

    ch_versions          = Channel.empty()

    //TODO: Temporary until the if's can be removed and printing to terminal is prevented with "when" in the modules.config
    deepvariant_vcf      = Channel.empty()
    freebayes_vcf        = Channel.empty()
    haplotypecaller_gvcf = Channel.empty()
    genotype_gvcf        = Channel.empty()
    manta_vcf            = Channel.empty()
    strelka_vcf          = Channel.empty()

    // Remap channel with intervals
    cram_recalibrated_intervals = cram_recalibrated.combine(intervals)
        .map{ meta, cram, crai, intervals ->
            sample = meta.sample
            new_intervals = intervals.baseName != "no_intervals" ? intervals : []
            id = new_intervals ? sample + "_" + new_intervals.baseName : sample
            [[ id: id, sample: meta.sample, gender: meta.gender, status: meta.status, patient: meta.patient ], cram, crai, new_intervals]
        }

    // Remap channel with gzipped intervals + indexes
    cram_recalibrated_intervals_gz_tbi = cram_recalibrated.combine(intervals_bed_gz_tbi)
        .map{ meta, cram, crai, bed, tbi ->
            sample = meta.sample
            new_bed = bed.simpleName != "no_intervals" ? bed : []
            new_tbi = tbi.simpleName != "no_intervals" ? tbi : []
            id = new_bed ? sample + "_" + new_bed.simpleName : sample
            new_meta = [ id: id, sample: meta.sample, gender: meta.gender, status: meta.status, patient: meta.patient ]
            [new_meta, cram, crai, new_bed, new_tbi]
        }

    // DEEPVARIANT
    if(params.tools.contains('deepvariant')){
        DEEPVARIANT(cram_recalibrated_intervals, fasta, fasta_fai, intervals_bed_combine_gz, num_intervals)

        deepvariant_vcf = DEEPVARIANT.out.deepvariant_vcf
        ch_versions     = ch_versions.mix(DEEPVARIANT.out.versions)
    }

    // FREEBAYES
    if (params.tools.contains('freebayes')){
        // Remap channel for Freebayes
        cram_recalibrated_intervals_freebayes = cram_recalibrated_intervals
            .map{ meta, cram, crai, intervals ->
                [meta, cram, crai, [], [], intervals]
            }
        FREEBAYES(cram_recalibrated_intervals_freebayes, fasta, fasta_fai)

        freebayes_vcf   = FREEBAYES.out.freebayes_vcf
        ch_versions     = ch_versions.mix(FREEBAYES.out.versions)
    }

    // HAPLOTYPECALLER
    if (params.tools.contains('haplotypecaller')){
        HAPLOTYPECALLER(cram_recalibrated_intervals,
                        fasta,
                        fasta_fai,
                        dict,
                        dbsnp,
                        dbsnp_tbi,
                        num_intervals,
                        intervals_bed_combine_gz,
                        intervals_bed_combine_gz_tbi)

        haplotypecaller_gvcf = HAPLOTYPECALLER.out.haplotypecaller_gvcf
        genotype_gvcf        = HAPLOTYPECALLER.out.genotype_gvcf
        ch_versions          = ch_versions.mix(HAPLOTYPECALLER.out.versions)

    }

    // MANTA
    if (params.tools.contains('manta')){
        MANTA_GERMLINE (cram_recalibrated_intervals_gz_tbi,
                        fasta,
                        fasta_fai,
                        intervals_bed_combine_gz,
                        num_intervals)

        manta_vcf   = MANTA_GERMLINE.out.manta_vcf
        ch_versions = ch_versions.mix(MANTA_GERMLINE.out.versions)
    }

    // STRELKA
    if (params.tools.contains('strelka')){
        STRELKA_SINGLE(cram_recalibrated_intervals_gz_tbi,
                fasta,
                fasta_fai,
                intervals_bed_combine_gz,
                num_intervals)

        strelka_vcf = STRELKA_SINGLE.out.strelka_vcf
        ch_versions = ch_versions.mix(STRELKA_SINGLE.out.versions)
    }

    //TIDDIT
    //TODO

    emit:
    deepvariant_vcf
    freebayes_vcf
    haplotypecaller_gvcf
    genotype_gvcf
    manta_vcf
    strelka_vcf

    versions = ch_versions
}
