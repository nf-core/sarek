//
// SMALL_GERMLINE_BENCHMARK: SUBWORKFLOW FOR SMALL GERMLINE VARIANTS
//

include { RTGTOOLS_VCFEVAL                                     } from '../../modules/nf-core/rtgtools/vcfeval/main'


workflow VCF_BENCHMARK_SMALL_VARIANTS {
    take:
    ch_test   // channel: test vcf coming from pipeline [val(meta), test.vcf.gz, test.vcf.gz.tbi]
    ch_truth  // channel: truth vcf [val(meta), truth.vcf.gz, truth.vcf.gz.tbi]
    ch_bed    // channel: bed file [val(meta), target.bed] //TODO: is optional for rtgvcfeval -> remove?

    main:
    versions        = Channel.empty()
    summary_reports = Channel.empty()

    // apply rtgtools eval method
    RTGTOOLS_VCFEVAL(
        // TODO: correct mapping according to input channels
        input_ch.map { meta, vcf, tbi, truth_vcf, truth_tbi, bed ->
            [ meta, vcf, tbi, truth_vcf, truth_tbi, bed, [] ]
        },
        [ [], [] ]
    )
    versions = versions.mix(RTGTOOLS_VCFEVAL.out.versions.first())

    // collect summary reports
    RTGTOOLS_VCFEVAL.out.summary
        .map { meta, file -> tuple([vartype: meta.vartype] + [benchmark_tool: "rtgtools"], file) }
        .groupTuple()
        .set{ report }

    summary_reports = summary_reports.mix(report)

    emit:
    versions
    summary_reports

}
