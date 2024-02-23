//
// CONTROLFREEC somatc variant calling
//
// For all modules here:
// A when clause condition is defined in the conf/modules.config to determine if the module should be run

include { CONTROLFREEC_FREEC              as FREEC_SOMATIC       } from '../../../modules/nf-core/controlfreec/freec/main'
include { CONTROLFREEC_ASSESSSIGNIFICANCE as ASSESS_SIGNIFICANCE } from '../../../modules/nf-core/controlfreec/assesssignificance/main'
include { CONTROLFREEC_FREEC2BED          as FREEC2BED           } from '../../../modules/nf-core/controlfreec/freec2bed/main'
include { CONTROLFREEC_FREEC2CIRCOS       as FREEC2CIRCOS        } from '../../../modules/nf-core/controlfreec/freec2circos/main'
include { CONTROLFREEC_MAKEGRAPH2         as MAKEGRAPH2           } from '../../../modules/nf-core/controlfreec/makegraph2/main'

workflow BAM_VARIANT_CALLING_SOMATIC_CONTROLFREEC {
    take:
    controlfreec_input       // channel: [mandatory] [meta, pileup_normal, pileup_tumor, [], [], [], []]
    fasta                    // channel: [mandatory]
    fasta_fai                // channel: [mandatory]
    dbsnp                    // channel: [mandatory]
    dbsnp_tbi                // channel: [mandatory]
    chr_files                // channel: [mandatory]
    mappability              // channel: [mandatory]
    intervals_bed            // channel: [optional]  Contains a bed file of all intervals combined provided with the cram input(s). Should be empty for WGS

    main:

    ch_versions = Channel.empty()

    FREEC_SOMATIC(controlfreec_input, fasta, fasta_fai, [], dbsnp, dbsnp_tbi, chr_files, mappability, intervals_bed, [])

    //Filter the files that come out of freec somatic as ASSESS_SIGNIFICANCE only takes one cnv and one ratio file
    //Creates empty channel if file is missing
    cnv_files = FREEC_SOMATIC.out.CNV
    .map{ meta, cnv ->
        def tumor_file = cnv instanceof List ? cnv.find { it.toString().endsWith("gz_CNVs") } : cnv //only find if its a list, else it returns only the filename without the path
        if (!tumor_file){
            log.error "CNVs tumor file not found for sample $meta.id"
        }
        [meta,tumor_file]
    }

    ratio_files = FREEC_SOMATIC.out.ratio
    .map{ meta, ratio ->
        def tumor_file = ratio instanceof List ? ratio.find { it.toString().endsWith("gz_ratio.txt") } : ratio //same here as cnv
        if (!tumor_file){
            log.error "Ratio tumor file not found for sample $meta.id"
        }
        [meta,tumor_file]
    }

    //Join the pairs
    assess_significance_input = cnv_files.join(ratio_files, failOnDuplicate: true, failOnMismatch: true)

    ASSESS_SIGNIFICANCE(assess_significance_input)
    FREEC2BED(FREEC_SOMATIC.out.ratio)
    FREEC2CIRCOS(FREEC_SOMATIC.out.ratio)
    MAKEGRAPH2(FREEC_SOMATIC.out.ratio.join(FREEC_SOMATIC.out.BAF, failOnDuplicate: true, failOnMismatch: true))

    ch_versions = ch_versions.mix(FREEC_SOMATIC.out.versions)
    ch_versions = ch_versions.mix(ASSESS_SIGNIFICANCE.out.versions)
    ch_versions = ch_versions.mix(FREEC2BED.out.versions)
    ch_versions = ch_versions.mix(FREEC2CIRCOS.out.versions)
    ch_versions = ch_versions.mix(MAKEGRAPH2.out.versions)

    emit:
    versions = ch_versions
}
