//
// GERMLINE VARIANT CALLING
//

include { RUN_DEEPVARIANT                           } from './variantcalling/deepvariant.nf'
include { RUN_FREEBAYES                             } from './variantcalling/freebayes.nf'
include { RUN_HAPLOTYPECALLER                       } from './variantcalling/haplotypecaller.nf'
include { RUN_MANTA                                 } from './variantcalling/manta.nf'
include { RUN_STRELKA                               } from './variantcalling/strelka.nf'
//include { RUN_TIDDIT                                } from './variantcalling/tiddit.nf'

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
        intervals_bed_combine_gz     // channel: [mandatory] intervals/target regions index zipped and indexed in one file
        num_intervals                // val: number of intervals that are used to parallelize exection, either based on capture kit or GATK recommended for WGS
        // joint_germline               // val: true/false on whether to run joint_germline calling, only works in combination with haplotypecaller at the moment

    main:

    ch_versions = Channel.empty()

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
    RUN_DEEPVARIANT(cram_recalibrated_intervals, fasta, fasta_fai, intervals_bed_combine_gz, num_intervals)

    // FREEBAYES
    // Remap channel for Freebayes
    cram_recalibrated_intervals_freebayes = cram_recalibrated_intervals
        .map{ meta, cram, crai, intervals ->
            [meta, cram, crai, [], [], intervals]
        }

    RUN_FREEBAYES(cram_recalibrated_intervals_freebayes, fasta, fasta_fai)

    // HAPLOTYPECALLER
    RUN_HAPLOTYPECALLER(cram_recalibrated_intervals,
                        fasta,
                        fasta_fai,
                        dict,
                        dbsnp,
                        dbsnp_tbi,
                        num_intervals)

    // MANTA
    RUN_MANTA(cram_recalibrated_intervals_gz_tbi,
        fasta,
        fasta_fai,
        num_intervals)

    // STRELKA
    RUN_STRELKA(cram_recalibrated_intervals_gz_tbi,
        fasta,
        fasta_fai,
        num_intervals)

    //TIDDIT
    //TODO

    ch_versions = ch_versions.mix(RUN_DEEPVARIANT.out.versions)
    ch_versions = ch_versions.mix(RUN_FREEBAYES.out.versions)
    ch_versions = ch_versions.mix(RUN_HAPLOTYPECALLER.out.versions)
    ch_versions = ch_versions.mix(RUN_MANTA.out.versions)
    ch_versions = ch_versions.mix(RUN_STRELKA.out.versions)

    emit:
    deepvariant_vcf = RUN_DEEPVARIANT.out.deepvariant_vcf
    freebayes_vcf   = RUN_FREEBAYES.out.freebayes_vcf
    haplotypecaller_gvcf = RUN_HAPLOTYPECALLER.out.haplotypecaller_gvcf
    genotype_gvcf   = RUN_HAPLOTYPECALLER.out.genotype_gvcf
    manta_vcf   = RUN_MANTA.out.manta_vcf
    strelka_vcf = RUN_STRELKA.out.strelka_vcf

    versions = ch_versions
}
