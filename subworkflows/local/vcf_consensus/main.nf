//
// Intersect VCFs and merge consensus variants from all callers
//

include { BCFTOOLS_ISEC        } from '../../../modules/nf-core/bcftools/isec'
include { BCFTOOLS_CONCAT      } from '../../../modules/nf-core/bcftools/concat'
include { CONSENSUS_FROM_SITES } from '../../../modules/local/consensus_from_sites'

workflow CONSENSUS {

    take:
    vcfs     // [meta, vcf ,tbi]

    main:
    ch_versions = Channel.empty()

    ch_vcfs = vcfs
        .branch{ meta, vcf, tbi ->
            // Somatic Strelka samples have tumor_id field (tumor-normal pairs)
            // This is semantically equivalent to checking status == '1' (tumor) but more explicit
            strelka_somatic: meta.variantcaller == 'strelka' && meta.tumor_id
            other: true
        }

    // Group somatic Strelka SNVs and INDELs by sample for concatenation
    // Remove filename from grouping key since SNVs and INDELs have different filenames but should be grouped together
    ch_strelka_grouped = ch_vcfs.strelka_somatic
        .map { meta, vcf, tbi ->
            def key = meta - meta.subMap('filename')
            [key, vcf, tbi]
        }
        .groupTuple(size: 2)

    BCFTOOLS_CONCAT(ch_strelka_grouped)// somatic strelkas have two vcf files: SNPs and indels
    ch_versions = ch_versions.mix(BCFTOOLS_CONCAT.out.versions)

    // Combine concat strelka with remaining VCFs
    // Bundle each VCF with its caller to preserve association through grouping
    ch_consensus_in = ch_vcfs.other
                        .mix(BCFTOOLS_CONCAT.out.vcf.join(BCFTOOLS_CONCAT.out.tbi))
                        .map { meta, vcf, tbi ->
                                    def caller = meta.variantcaller
                                    def groupKey = meta - meta.subMap('variantcaller', 'contamination', 'filename', 'data_type', 'num_intervals')
                                    [groupKey, [vcf, caller], tbi]
                        }
                        // TODO: blocking operation unless we learn how many variantcallers were
                        // specified - also depends on whether this is n, t, or nt and how many
                        // variantcallers are actually executed
                        .groupTuple()
                        .map { meta, vcf_caller_pairs, tbis ->
                            // Sort by vcf name for predictable isec input order
                            // callers list will match isec output order in sites.txt
                            def sorted_pairs = vcf_caller_pairs.sort { a, b -> a[0].name <=> b[0].name }
                            def sorted_vcfs = sorted_pairs.collect { it[0] }
                            def callers = sorted_pairs.collect { it[1] }
                            [meta + [callers: callers], sorted_vcfs, tbis]
                        }


    BCFTOOLS_ISEC(ch_consensus_in)
    ch_versions = ch_versions.mix(BCFTOOLS_ISEC.out.versions)

    // Filter out empty isec results (no consensus variants found)
    ch_isec_with_results = BCFTOOLS_ISEC.out.results
        .filter { meta, dir ->
            def sites_file = dir.resolve('sites.txt')
            sites_file.exists() && sites_file.size() > 0
        }

    // Create consensus VCF from sites.txt with caller presence info
    // Versions are collected via topic channel
    CONSENSUS_FROM_SITES(ch_isec_with_results)

    emit:
    versions = ch_versions
    vcfs     = CONSENSUS_FROM_SITES.out.vcf
    tbis     = CONSENSUS_FROM_SITES.out.tbi

}
