//
// TUMOR VARIANT CALLING
// Should be only run on patients without normal sample
//

//include { RUN_CONTROLFREEC                        } from '../nf-core/variantcalling/controlfreec/main.nf'
include { RUN_FREEBAYES                           } from '../nf-core/variantcalling/freebayes/main.nf'
include { GATK_TUMOR_ONLY_SOMATIC_VARIANT_CALLING } from '../../subworkflows/nf-core/gatk4/tumor_only_somatic_variant_calling/main'
include { RUN_MANTA_TUMORONLY                     } from '../nf-core/variantcalling/manta/tumoronly/main.nf'
include { RUN_STRELKA_SINGLE                      } from '../nf-core/variantcalling/strelka/single/main.nf'
include { RUN_CONTROLFREEC                        } from '../nf-core/variantcalling/controlfreec/somatic/main.nf'

workflow TUMOR_ONLY_VARIANT_CALLING {
    take:
        tools                        // Mandatory, list of tools to apply
        cram_recalibrated            // channel: [mandatory] cram
        dbsnp                        // channel: [mandatory] dbsnp
        dbsnp_tbi                    // channel: [mandatory] dbsnp_tbi
        dict                         // channel: [mandatory] dict
        fasta                        // channel: [mandatory] fasta
        fasta_fai                    // channel: [mandatory] fasta_fai
        intervals                    // channel: [mandatory] intervals/target regions
        intervals_bed_gz_tbi         // channel: [mandatory] intervals/target regions index zipped and indexed
        intervals_bed_combine_gz_tbi // channel: [mandatory] intervals/target regions index zipped and indexed
        intervals_bed_combine_gz     // channel: [mandatory] intervals/target regions index zipped and indexed in one file
        intervals_bed_combined        // channel: [mandatory] intervals/target regions in one file unzipped
        num_intervals                // val: number of intervals that are used to parallelize exection, either based on capture kit or GATK recommended for WGS
        no_intervals
        germline_resource            // channel: [optional]  germline_resource
        germline_resource_tbi        // channel: [optional]  germline_resource_tbi
        panel_of_normals             // channel: [optional]  panel_of_normals
        panel_of_normals_tbi         // channel: [optional]  panel_of_normals_tbi
        mappability

    main:

    ch_versions         = Channel.empty()

    //TODO: Temporary until the if's can be removed and printing to terminal is prevented with "when" in the modules.config
    freebayes_vcf       = Channel.empty()
    manta_vcf           = Channel.empty()
    mutect2_vcf         = Channel.empty()
    strelka_vcf         = Channel.empty()

    cram_recalibrated.combine(intervals).map{ meta, cram, crai, intervals ->
        sample = meta.sample
        new_intervals = intervals.baseName != "no_intervals" ? intervals : []
        id = new_intervals ? sample + "_" + new_intervals.baseName : sample
        new_new_meta = [ id: id, sample: meta.sample, gender: meta.gender, status: meta.status, patient: meta.patient ]
        [new_new_meta, cram, crai, new_intervals]
    }.set{cram_recalibrated_intervals}

    cram_recalibrated.combine(intervals_bed_gz_tbi)
        .map{ meta, cram, crai, bed, tbi ->
            sample = meta.sample
            new_bed = bed.simpleName != "no_intervals" ? bed : []
            new_tbi = tbi.simpleName != "no_intervals" ? tbi : []
            id = new_bed ? sample + "_" + new_bed.simpleName : sample
            new_meta = [ id: id, sample: meta.sample, gender: meta.gender, status: meta.status, patient: meta.patient ]
            [new_meta, cram, crai, new_bed, new_tbi]
        }.set{cram_recalibrated_intervals_gz_tbi}

    if(tools.contains('controlfreec')){
        cram_recalibrated_intervals.map {meta, cram, crai, intervals -> [meta, cram, intervals]}.set{cram_intervals_no_index}
        RUN_CONTROLFREEC(cram_intervals_no_index,
                        fasta,
                        fasta_fai,
                        dbsnp,
                        dbsnp_tbi,
                        chr_files,
                        mappability,
                        intervals_bed_combined,
                        num_intervals)
        ch_versions = ch_versions.mix(RUN_CONTROLFREEC.out.versions)
    }

    if (tools.contains('freebayes')){
        // Remap channel for Freebayes
        cram_recalibrated_intervals_freebayes = cram_recalibrated_intervals
            .map{ meta, cram, crai, intervals ->
                [meta, cram, crai, [], [], intervals]
            }

        RUN_FREEBAYES(cram_recalibrated_intervals_freebayes, fasta, fasta_fai, intervals_bed_combine_gz, num_intervals)

        freebayes_vcf = RUN_FREEBAYES.out.freebayes_vcf
        ch_versions   = ch_versions.mix(RUN_FREEBAYES.out.versions)
    }

    if (tools.contains('mutect2')) {

        which_norm = []
        cram_recalibrated_intervals.map{ meta, cram, crai, intervals -> [meta, cram, crai, intervals, which_norm]}.set{cram_recalibrated_mutect2}
        GATK_TUMOR_ONLY_SOMATIC_VARIANT_CALLING(cram_recalibrated_mutect2,
                                                fasta,
                                                fasta_fai,
                                                dict,
                                                germline_resource,
                                                germline_resource_tbi,
                                                panel_of_normals,
                                                panel_of_normals_tbi,
                                                intervals_bed_combine_gz,
                                                num_intervals)

        mutect2_vcf = GATK_TUMOR_ONLY_SOMATIC_VARIANT_CALLING.out.mutect2_vcf
        ch_versions = ch_versions.mix(GATK_TUMOR_ONLY_SOMATIC_VARIANT_CALLING.out.versions)

    }

    if (tools.contains('manta')){
        //TODO: Research if splitting by intervals is ok, we pretend for now it is fine. Seems to be the consensus on upstream modules implementaiton too

        RUN_MANTA_TUMORONLY(cram_recalibrated_intervals_gz_tbi,
                            fasta,
                            fasta_fai,
                            intervals_bed_combine_gz,
                            num_intervals)

        manta_vcf   = RUN_MANTA_TUMORONLY.out.manta_vcf
        ch_versions = ch_versions.mix(RUN_MANTA_TUMORONLY.out.versions)
    }

    if (tools.contains('strelka')) {
        RUN_STRELKA_SINGLE( cram_recalibrated_intervals_gz_tbi,
                            fasta,
                            fasta_fai,
                            intervals_bed_combine_gz,
                            num_intervals)

        strelka_vcf = RUN_STRELKA_SINGLE.out.strelka_vcf
        ch_versions = ch_versions.mix(RUN_STRELKA_SINGLE.out.versions)
    }


    // if (tools.contains('tiddit')){
    // }

    emit:
    versions = ch_versions

    freebayes_vcf
    manta_vcf
    mutect2_vcf
    strelka_vcf

}
