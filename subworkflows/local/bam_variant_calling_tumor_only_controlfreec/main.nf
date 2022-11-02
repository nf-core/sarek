include { CONTROLFREEC_FREEC as FREEC_TUMORONLY                  } from '../../../modules/nf-core/controlfreec/freec/main'
include { CONTROLFREEC_ASSESSSIGNIFICANCE as ASSESS_SIGNIFICANCE } from '../../../modules/nf-core/controlfreec/assesssignificance/main'
include { CONTROLFREEC_FREEC2BED as FREEC2BED                    } from '../../../modules/nf-core/controlfreec/freec2bed/main'
include { CONTROLFREEC_FREEC2CIRCOS as FREEC2CIRCOS              } from '../../../modules/nf-core/controlfreec/freec2circos/main'
include { CONTROLFREEC_MAKEGRAPH as MAKEGRAPH                    } from '../../../modules/nf-core/controlfreec/makegraph/main'

workflow BAM_VARIANT_CALLING_TUMOR_ONLY_CONTROLFREEC {
    take:
    controlfreec_input       // channel: [mandatory] [meta, [], pileup_tumor, [], [], [], []]
    fasta                    // channel: [mandatory]
    fasta_fai                // channel: [mandatory]
    dbsnp                    // channel: [mandatory]
    dbsnp_tbi                // channel: [mandatory]
    chr_files                // channel: [mandatory]
    mappability              // channel: [mandatory]
    intervals_bed            // channel: [optional]  Contains a bed file of all intervals combined provided with the cram input(s). Should be empty for WGS

    main:

    ch_versions = Channel.empty()

    FREEC_TUMORONLY(controlfreec_input,
                fasta,
                fasta_fai,
                [],
                dbsnp,
                dbsnp_tbi,
                chr_files,
                mappability,
                intervals_bed,
                [])

    ASSESS_SIGNIFICANCE(FREEC_TUMORONLY.out.CNV.join(FREEC_TUMORONLY.out.ratio))
    FREEC2BED(FREEC_TUMORONLY.out.ratio)
    FREEC2CIRCOS(FREEC_TUMORONLY.out.ratio)
    MAKEGRAPH(FREEC_TUMORONLY.out.ratio.join(FREEC_TUMORONLY.out.BAF))

    ch_versions = ch_versions.mix(FREEC_TUMORONLY.out.versions)
    ch_versions = ch_versions.mix(ASSESS_SIGNIFICANCE.out.versions)
    ch_versions = ch_versions.mix(FREEC2BED.out.versions)
    ch_versions = ch_versions.mix(FREEC2CIRCOS.out.versions)
    ch_versions = ch_versions.mix(MAKEGRAPH.out.versions)

    emit:
    versions = ch_versions
}
