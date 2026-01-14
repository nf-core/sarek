//
// Intersect VCFs
//

include { BCFTOOLS_ISEC } from '../../../modules/nf-core/bcftools/isec'
include { BCFTOOLS_CONCAT } from '../../../modules/nf-core/bcftools/concat'

workflow CONSENSUS {

    take:
    vcfs     // [meta, vcf ,tbi]

    main:
    ch_versions = Channel.empty()

    ch_vcfs = vcfs
        .branch{ meta, vcf, tbi ->
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

    //Combine concat strelka with remaining VCFs
    ch_consensus_in = ch_vcfs.other
                        .mix(BCFTOOLS_CONCAT.out.vcf.join(BCFTOOLS_CONCAT.out.tbi))
                        .map { meta, vcf, tbi ->
                                    [meta - meta.subMap('variantcaller', 'contamination', 'filename'), vcf, tbi]
                        }
                        //TODO blocking operation unless we learn how many variantcallers were
                        // specified also this depends on whether this n,t, or nt on how many
                        //variantcallers are actually executed
                        .groupTuple()
                        .map {meta, vcf, tbi ->
                            // Sorting the VCF files to ensure the consensus calling is done in a predictable manner
                            def vcf_sorted = (vcf instanceof List) ? vcf.sort() : vcf
                            [meta, vcf_sorted, tbi]
                        }


    BCFTOOLS_ISEC(ch_consensus_in)
    ch_versions = ch_versions.mix(BCFTOOLS_ISEC.out.versions)

    ch_consensus_results = BCFTOOLS_ISEC.out.results
        .map { meta, dir ->
            def files = dir.listFiles()
            def vcf = files.find { it.name == '0000.vcf.gz' }
            def tbi = files.find { it.name == '0000.vcf.gz.tbi' }
            // return the consensus file for annotation
            return [meta, vcf, tbi]
        }

    emit:
    versions = ch_versions
    vcfs = ch_consensus_results.map{ meta, vcf_, _tbi -> [meta, vcf_]}
    tbis = ch_consensus_results.map{ meta, _vcf, tbi_ -> [meta, tbi_]}

}
