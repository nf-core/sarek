//
// GERMLINE VARIANT CALLING
//

include { RUN_DEEPVARIANT     } from '../nf-core/variantcalling/deepvariant/main.nf'
include { RUN_FREEBAYES       } from '../nf-core/variantcalling/freebayes/main.nf'
include { RUN_HAPLOTYPECALLER } from '../nf-core/variantcalling/haplotypecaller/main.nf'
include { RUN_MANTA_GERMLINE  } from '../nf-core/variantcalling/manta/germline/main.nf'
include { RUN_STRELKA_SINGLE  } from '../nf-core/variantcalling/strelka/single/main.nf'
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
    haplotypecaller_vcf  = Channel.empty()
    genotype_gvcf        = Channel.empty()
    manta_vcf            = Channel.empty()
    strelka_vcf          = Channel.empty()

    // Remap channel with intervals
    cram_recalibrated_intervals = cram_recalibrated.combine(intervals)
        .map{ meta, cram, crai, intervals ->
            sample = meta.sample
            //new_intervals = num_intervals > 1 ? intervals : []
            new_intervals = intervals.baseName != "no_intervals" ? intervals : []
            id = new_intervals ? sample + "_" + new_intervals.baseName : sample
            [[ id: id, sample: meta.sample, gender: meta.gender, status: meta.status, patient: meta.patient ], cram, crai, new_intervals]
        }

    // Remap channel with gzipped intervals + indexes
    cram_recalibrated_intervals_gz_tbi = cram_recalibrated.combine(intervals_bed_gz_tbi)
        .map{ meta, cram, crai, bed, tbi ->
            sample = meta.sample
            //new_bed = num_intervals > 1 ? bed : [] //TODO can I pass in empty lists? Then I only need to work with the id line
            new_bed = bed.simpleName != "no_intervals" ? bed : []
            new_tbi = tbi.simpleName != "no_intervals" ? tbi : []
            id = new_bed ? sample + "_" + new_bed.simpleName : sample
            new_meta = [ id: id, sample: meta.sample, gender: meta.gender, status: meta.status, patient: meta.patient ]
            [new_meta, cram, crai, new_bed, new_tbi]
        }

    // DEEPVARIANT
    if(params.tools.contains('deepvariant')){
        RUN_DEEPVARIANT(cram_recalibrated_intervals, fasta, fasta_fai, intervals_bed_combine_gz, num_intervals)

        deepvariant_vcf = RUN_DEEPVARIANT.out.deepvariant_vcf
        ch_versions     = ch_versions.mix(RUN_DEEPVARIANT.out.versions)
    }

    // FREEBAYES
    if (params.tools.contains('freebayes')){
        // Remap channel for Freebayes
        cram_recalibrated_intervals_freebayes = cram_recalibrated_intervals
            .map{ meta, cram, crai, intervals ->
                [meta, cram, crai, [], [], intervals]
            }
        RUN_FREEBAYES(cram_recalibrated_intervals_freebayes, fasta, fasta_fai, intervals_bed_combine_gz, num_intervals)

        freebayes_vcf   = RUN_FREEBAYES.out.freebayes_vcf
        ch_versions     = ch_versions.mix(RUN_FREEBAYES.out.versions)
    }

    // HAPLOTYPECALLER
    if (params.tools.contains('haplotypecaller')){
        RUN_HAPLOTYPECALLER(cram_recalibrated_intervals,
                        fasta,
                        fasta_fai,
                        dict,
                        dbsnp,
                        dbsnp_tbi,
                        intervals_bed_combine_gz,
                        intervals_bed_combine_gz_tbi,
                        num_intervals)

        haplotypecaller_vcf  = RUN_HAPLOTYPECALLER.out.haplotypecaller_vcf
        //genotype_gvcf        = RUN_HAPLOTYPECALLER.out.genotype_gvcf
        ch_versions          = ch_versions.mix(RUN_HAPLOTYPECALLER.out.versions)
    }

    // MANTA
    if (params.tools.contains('manta')){
        RUN_MANTA_GERMLINE (cram_recalibrated_intervals_gz_tbi,
                        fasta,
                        fasta_fai,
                        intervals_bed_combine_gz,
                        num_intervals)

        manta_vcf   = RUN_MANTA_GERMLINE.out.manta_vcf
        ch_versions = ch_versions.mix(RUN_MANTA_GERMLINE.out.versions)
    }

    // STRELKA
    if (params.tools.contains('strelka')){
        RUN_STRELKA_SINGLE(cram_recalibrated_intervals_gz_tbi,
                fasta,
                fasta_fai,
                intervals_bed_combine_gz,
                num_intervals)

        strelka_vcf = RUN_STRELKA_SINGLE.out.strelka_vcf
        ch_versions = ch_versions.mix(RUN_STRELKA_SINGLE.out.versions)
    }

    //TIDDIT
    //TODO

    emit:
    deepvariant_vcf
    freebayes_vcf
    haplotypecaller_vcf
    genotype_gvcf
    manta_vcf
    strelka_vcf

    versions = ch_versions
}
