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
        .map{ meta, cram, crai, intervals, num_intervals ->
            // If either no scatter/gather is done, i.e. no interval (0) or one interval (1), then don't rename samples
            def new_id = num_intervals <= 1 ? meta.sample : meta.sample + "_" + intervals.baseName

            //If no interval file provided (0) then add empty list
            intervals_new = num_intervals == 0 ? [] : intervals

            [[patient:meta.patient, sample:meta.sample, gender:meta.gender, status:meta.status, id:new_id, data_type:meta.data_type, num_intervals:num_intervals],
            cram, crai, intervals_new]
        }

    // Remap channel with gzipped intervals + indexes
    cram_recalibrated_intervals_gz_tbi = cram_recalibrated.combine(intervals_bed_gz_tbi)
        .map{ meta, cram, crai, bed_tbi, num_intervals ->
            // If either no scatter/gather is done, i.e. no interval (0) or one interval (1), then don't rename samples
            def new_id = num_intervals <= 1 ? meta.sample : meta.sample + "_" + bed_tbi[0].simpleName

            //If no interval file provided (0) then add empty list
            bed_new = num_intervals == 0 ? [] : bed_tbi[0]
            tbi_new = num_intervals == 0 ? [] : bed_tbi[1]

            [[patient:meta.patient, sample:meta.sample, gender:meta.gender, status:meta.status, id:new_id, data_type:meta.data_type, num_intervals:num_intervals],
            cram, crai, bed_new, tbi_new]
        }

    // DEEPVARIANT
    if(params.tools.contains('deepvariant')){
        RUN_DEEPVARIANT(cram_recalibrated_intervals, fasta, fasta_fai, intervals_bed_combine_gz)

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
        RUN_FREEBAYES(cram_recalibrated_intervals_freebayes, fasta, fasta_fai, intervals_bed_combine_gz)

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
                        intervals_bed_combine_gz_tbi)

        haplotypecaller_vcf  = RUN_HAPLOTYPECALLER.out.haplotypecaller_vcf
        //genotype_gvcf        = RUN_HAPLOTYPECALLER.out.genotype_gvcf
        ch_versions          = ch_versions.mix(RUN_HAPLOTYPECALLER.out.versions)
    }

    // MANTA
    if (params.tools.contains('manta')){
        RUN_MANTA_GERMLINE (cram_recalibrated_intervals_gz_tbi,
                        fasta,
                        fasta_fai,
                        intervals_bed_combine_gz)

        manta_vcf   = RUN_MANTA_GERMLINE.out.manta_vcf
        ch_versions = ch_versions.mix(RUN_MANTA_GERMLINE.out.versions)
    }

    // STRELKA
    if (params.tools.contains('strelka')){
        RUN_STRELKA_SINGLE(cram_recalibrated_intervals_gz_tbi,
                dict,
                fasta,
                fasta_fai,
                intervals_bed_combine_gz)

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
