//
// RUN BENCHMARKING OF SMALL GERMLINE VARIANTS
// Taken from nf-core/variantbenchmarking
//

include { RTGTOOLS_FORMAT               } from '../../../modules/nf-core/rtgtools/format/main'
include { RTGTOOLS_VCFEVAL              } from '../../../modules/nf-core/rtgtools/vcfeval/main'
include { BCFTOOLS_INDEX as INDEX_QUERY } from '../../../modules/nf-core/bcftools/index/main.nf'
include { BCFTOOLS_INDEX as INDEX_TRUTH } from '../../../modules/nf-core/bcftools/index/main.nf'

workflow VCF_BENCHMARK {

    take:
    input_ch  // channel: [val(meta),test_vcf]
    input_truth // channel: [val(meta),truth_vcf,truth_bed]
    fasta     // reference channel [val(meta), ref.fa]
    fai       // reference channel [val(meta), ref.fa.fai]

    main:
    versions        = Channel.empty()
    summary_reports = Channel.empty()

    query_vcf_tbi = Channel.empty()
    truth_vcf_tbi = Channel.empty()
    INDEX_TRUTH(input_ch.map { meta, test_vcf, truth_vcf, bed -> [ meta, test_vcf  ] })
    INDEX_QUERY(input_ch.map { meta, test_vcf, truth_vcf, bed -> [ meta, truth_vcf ] })

    versions = versions.mix(INDEX_TRUTH.out.versions)
    versions = versions.mix(INDEX_QUERY.out.versions)

    query_vcf_tbi = query_vcf_tbi.mix(INDEX_QUERY.out.tbi)
    truth_vcf_tbi = truth_vcf_tbi.mix(INDEX_TRUTH.out.tbi)

    // Use rtgtools format to generate sdf file if necessary
    RTGTOOLS_FORMAT(
        fasta.map { meta, fasta -> [ meta, fasta, [], [] ] }
    )
    versions = versions.mix(RTGTOOLS_FORMAT.out.versions)
    sdf = RTGTOOLS_FORMAT.out.sdf

    // apply rtgtools eval method
    RTGTOOLS_VCFEVAL(
        input_ch.map { meta, query_vcf, query_vcf_tbi, truth_vcf, truth_vcf_tbi, bed ->
            [ meta, query_vcf, query_vcf_tbi, truth_vcf, truth_vcf_tbi, bed, [] ]
        },
        sdf
    )
    versions = versions.mix(RTGTOOLS_VCFEVAL.out.versions.first())

    // collect summary reports
    RTGTOOLS_VCFEVAL.out.summary
        .map { meta, file -> tuple([vartype: meta.vartype] + [benchmark_tool: "rtgtools"], file) }
        .groupTuple()
        .set{ report }

    summary_reports = summary_reports.mix(report)

    emit:
    summary_reports
    versions

    }
