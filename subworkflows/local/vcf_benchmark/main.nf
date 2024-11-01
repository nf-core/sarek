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
    input_truth // channel: [val(meta),truth_vcf]
    truth_bed // channel: [val(meta),truth_bed]
    fasta     // reference channel [val(meta), ref.fa]
    fai       // reference channel [val(meta), ref.fa.fai]

    main:
    input_ch.view()
    versions        = Channel.empty()
    summary_reports = Channel.empty()

    query_vcf_tbi = Channel.empty()
    truth_vcf_tbi = Channel.empty()

    INDEX_TRUTH( input_truth )
    INDEX_QUERY( input_ch    )

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

   // Combine input_ch and input_truth with query_vcf_tbi and truth_vcf_tbi
    combined_inputs = input_ch.combine(input_truth)
        .combine(truth_bed)
        .combine(query_vcf_tbi)
        .combine(truth_vcf_tbi)

    // Apply rtgtools eval method
    RTGTOOLS_VCFEVAL(
        combined_inputs.map { ch, truth, bed, query_tbi, truth_tbi ->
        def (meta, query_vcf) = ch
        def (meta_truth, truth_vcf) = truth
        def (meta_bed, truth_bed) = bed
        [ meta, query_vcf, query_tbi, truth_vcf, truth_tbi, truth_bed ]
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
