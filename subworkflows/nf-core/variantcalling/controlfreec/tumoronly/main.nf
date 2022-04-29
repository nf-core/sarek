include { CAT_CAT as CAT_MPILEUP_TUMOR                           } from '../../../../../modules/nf-core/modules/cat/cat/main.nf'
include { CONTROLFREEC_FREEC as FREEC_TUMORONLY                  } from '../../../../../modules/nf-core/modules/controlfreec/freec/main'
include { CONTROLFREEC_ASSESSSIGNIFICANCE as ASSESS_SIGNIFICANCE } from '../../../../../modules/nf-core/modules/controlfreec/assesssignificance/main'
include { CONTROLFREEC_FREEC2BED as FREEC2BED                    } from '../../../../../modules/nf-core/modules/controlfreec/freec2bed/main'
include { CONTROLFREEC_FREEC2CIRCOS as FREEC2CIRCOS              } from '../../../../../modules/nf-core/modules/controlfreec/freec2circos/main'
include { CONTROLFREEC_MAKEGRAPH as MAKEGRAPH                    } from '../../../../../modules/nf-core/modules/controlfreec/makegraph/main'
include { SAMTOOLS_MPILEUP as MPILEUP_TUMOR                      } from '../../../../../modules/nf-core/modules/samtools/mpileup/main'

workflow RUN_CONTROLFREEC_TUMORONLY {
    take:
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

    MPILEUP_TUMOR(cram_tumor, fasta)

    MPILEUP_TUMOR.out.mpileup.branch{
            intervals:    num_intervals > 1
            no_intervals: num_intervals == 1
        }.set{mpileup_tumor}

    //Merge mpileup only when intervals and natural order sort them
    CAT_MPILEUP_TUMOR(mpileup_tumor.intervals
        .map{ meta, pileup ->
            new_meta = meta.clone()
            new_meta.id = new_meta.sample
            [new_meta, pileup]
        }
        .groupTuple(size: num_intervals, sort:true))

    controlfreec_input_tumor = Channel.empty().mix(
        CAT_MPILEUP_TUMOR.out.file_out,
        mpileup_tumor.no_intervals
    ).map{ meta, pileup ->
        new_meta = meta.clone()
        new_meta.id = new_meta.sample
        [new_meta, pileup]
    }

    controlfreec_input_tumor
        .map{ meta, pileup_tumor ->
            [meta, [], pileup_tumor, [], [], [], []]
        }.set{controlfreec_input}

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

    ASSESS_SIGNIFICANCE( FREEC_TUMORONLY.out.CNV.join(FREEC_TUMORONLY.out.ratio))
    FREEC2BED( FREEC_TUMORONLY.out.ratio )
    FREEC2CIRCOS( FREEC_TUMORONLY.out.ratio )
    MAKEGRAPH(FREEC_TUMORONLY.out.ratio.join(FREEC_TUMORONLY.out.BAF))

    ch_versions = ch_versions.mix(MPILEUP_TUMOR.out.versions)
    ch_versions = ch_versions.mix(CAT_MPILEUP_TUMOR.out.versions)
    ch_versions = ch_versions.mix(FREEC_TUMORONLY.out.versions)
    ch_versions = ch_versions.mix(ASSESS_SIGNIFICANCE.out.versions)
    ch_versions = ch_versions.mix(FREEC2BED.out.versions)
    ch_versions = ch_versions.mix(FREEC2CIRCOS.out.versions)
    ch_versions = ch_versions.mix(MAKEGRAPH.out.versions)

    emit:
    versions = ch_versions
}
