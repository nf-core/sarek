include { CAT_CAT as CAT_MPILEUP_NORMAL                          } from '../../../../../modules/nf-core/modules/cat/cat/main.nf'
include { CAT_CAT as CAT_MPILEUP_TUMOR                           } from '../../../../../modules/nf-core/modules/cat/cat/main.nf'
include { CONTROLFREEC_FREEC as FREEC_SOMATIC                    } from '../../../../../modules/nf-core/modules/controlfreec/freec/main'
include { CONTROLFREEC_ASSESSSIGNIFICANCE as ASSESS_SIGNIFICANCE } from '../../../../../modules/nf-core/modules/controlfreec/assesssignificance/main'
include { CONTROLFREEC_FREEC2BED as FREEC2BED                    } from '../../../../../modules/nf-core/modules/controlfreec/freec2bed/main'
include { CONTROLFREEC_FREEC2CIRCOS as FREEC2CIRCOS              } from '../../../../../modules/nf-core/modules/controlfreec/freec2circos/main'
include { CONTROLFREEC_MAKEGRAPH as MAKEGRAPH                    } from '../../../../../modules/nf-core/modules/controlfreec/makegraph/main'
include { MPILEUP                                                } from '../../../../../subworkflows/nf-core/variantcalling/mpileup/main'

workflow RUN_CONTROLFREEC_SOMATIC {
    take:
    cram_normal              // channel: [mandatory] [meta, cram, crai, interval]
    cram_tumor               // channel: [mandatory] [meta, cram, crai, interval]
    fasta                    // channel: [mandatory]
    fasta_fai                // channel: [mandatory]
    dbsnp                    // channel: [mandatory]
    dbsnp_tbi                // channel: [mandatory]
    chr_files                // channel: [mandatory]
    mappability              // channel: [mandatory]
    intervals_bed            // channel: [optional]  Contains a bed file of all intervals combined provided with the cram input(s). Should be empty for WGS

    main:

    ch_versions = Channel.empty()

    MPILEUP(cram_normal, cram_tumor, fasta)

    controlfreec_input_normal = Channel.empty().mix(
        MPILEUP.out.cat_mpileup_normal,
        MPILEUP.out.mpileup_normal_no_intervals
    ).map{ meta, pileup ->
        new_meta = [patient:meta.patient, normal_id:meta.normal_id, tumor_id:meta.tumor_id, gender:meta.gender, id:meta.tumor_id + "_vs_" + meta.normal_id, num_intervals:meta.num_intervals]

        [new_meta, pileup]
    }

    controlfreec_input_tumor = Channel.empty().mix(
        MPILEUP.out.cat_mpileup_tumor,
        MPILEUP.out.mpileup_tumor_no_intervals
    ).map{ meta, pileup ->
        new_meta = [patient:meta.patient, normal_id:meta.normal_id, tumor_id:meta.tumor_id, gender:meta.gender, id:meta.tumor_id + "_vs_" + meta.normal_id, num_intervals:meta.num_intervals]
        [new_meta, pileup]
    }

    controlfreec_input_normal.cross(controlfreec_input_tumor)
        .map{ normal, tumor ->
            [normal[0], normal[1], tumor[1], [], [], [], []]
        }
        .set{controlfreec_input}

    FREEC_SOMATIC(controlfreec_input,
                fasta,
                fasta_fai,
                [],
                dbsnp,
                dbsnp_tbi,
                chr_files,
                mappability,
                intervals_bed,
                [])

    ASSESS_SIGNIFICANCE( FREEC_SOMATIC.out.CNV.join(FREEC_SOMATIC.out.ratio))
    FREEC2BED( FREEC_SOMATIC.out.ratio )
    FREEC2CIRCOS( FREEC_SOMATIC.out.ratio )
    MAKEGRAPH(FREEC_SOMATIC.out.ratio.join(FREEC_SOMATIC.out.BAF))

    ch_versions = ch_versions.mix(MPILEUP.out.versions)
    ch_versions = ch_versions.mix(FREEC_SOMATIC.out.versions)
    ch_versions = ch_versions.mix(ASSESS_SIGNIFICANCE.out.versions)
    ch_versions = ch_versions.mix(FREEC2BED.out.versions)
    ch_versions = ch_versions.mix(FREEC2CIRCOS.out.versions)
    ch_versions = ch_versions.mix(MAKEGRAPH.out.versions)

    emit:
    versions = ch_versions
}
