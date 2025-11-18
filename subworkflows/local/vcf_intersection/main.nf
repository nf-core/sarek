//
// Intersect VCFs
//

include { BCFTOOLS_ISEC } from '../../../modules/nf-core/bcftools/isec'
include { BCFTOOLS_CONCAT } from '../../../modules/nf-core/bcftools/concat'

workflow INTERSECTION {

    take:
    vcfs     // [meta, vcf ,tbi]

    main:
    ch_versions = Channel.empty()

    ch_vcfs = vcfs
        .branch{ meta, vcf, tbi ->
            strelka_somatic: meta.variantcaller == 'strelka' && meta.status == '1'
            other: true
        }

    BCFTOOLS_CONCAT(ch_vcfs.strelka_somatic.groupTuple(size: 2))// somatic strelkas have two vcf files: SNPs and indels
    ch_versions = ch_versions.mix(BCFTOOLS_CONCAT.out.versions)

    //Combine concat strelka with remaining VCFs
    ch_intersect_in = ch_vcfs.other
                        .mix(BCFTOOLS_CONCAT.out.vcf.join(BCFTOOLS_CONCAT.out.tbi))
                        .map { meta, vcf, tbi ->
                                    [meta - meta.subMap('variantcaller', 'contamination', 'filename'), vcf, tbi]
                        }
                        //TODO blocking operation unless we learn how many variantcallers were
                        // specified also this depends on whether this n,t, or nt on how many
                        //variantcallers are actually executed
                        .groupTuple()
                        .map {meta, vcf, tbi ->
                            // Sorting the VCF files to ensure the intersection is done in a predictable manner
                            def vcf_sorted = (vcf instanceof List) ? vcf.sort() : vcf
                            [meta, vcf_sorted, tbi]
                        }


    BCFTOOLS_ISEC(ch_intersect_in)
    ch_versions = ch_versions.mix(BCFTOOLS_ISEC.out.versions)

    ch_intersect_results = BCFTOOLS_ISEC.out.results
        .map { meta, dir ->
            def files = dir.listFiles()
            def vcf = files.find { it.name == '0000.vcf.gz' }
            def tbi = files.find { it.name == '0000.vcf.gz.tbi' }
            // return the intersected file for annotation
            return [meta, vcf, tbi]
        }

    emit:
    versions = ch_versions
    vcfs = ch_intersect_results.map{ meta, vcf_, _tbi -> [meta, vcf_]}
    tbis = ch_intersect_results.map{ meta, _vcf, tbi_ -> [meta, tbi_]}

}
