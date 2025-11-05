//
// Intersect VCFs
//

include { BCFTOOLS_ISEC } from '../../../modules/nf-core/bcftools/isec'
include { BCFTOOLS_CONCAT } from '../../../modules/nf-core/bcftools/concat'
workflow INTERSECTION {

    take:
    vcfs     // TODO ensure the input vcfs are sorted, should be of format vcf ,tbi

    main:

    versions = Channel.empty()

    vcfs.view()

    ch_vcfs = vcfs
        .branch{ meta, vcf, tbi ->
            strelka_somatic: meta.variantcaller == 'strelka' && meta.status == '1'
            other: true
        }

    BCFTOOLS_CONCAT(ch_vcfs.strelka_somatic.groupTuple(size: 2))// somatic strelkas have two vcf files: SNPs and indels
    versions = versions.mix(BCFTOOLS_CONCAT.out.versions)

    //Combine concat strelka with remaining VCFs
    ch_intersect_in = ch_vcfs.other
                        .mix(BCFTOOLS_CONCAT.out.vcf.join(BCFTOOLS_CONCAT.out.tbi))
                        .map { meta, vcf, tbi ->
                                    // TODO think about how much we want to keep here, dropping entire meta map right now
                                    [[id:meta.id], vcf, tbi]
                                }.groupTuple() //blocking operation unless we learn how many variantcallers were specified

    BCFTOOLS_ISEC(ch_intersect_in)
    versions = versions.mix(BCFTOOLS_ISEC.out.versions)

    //TODO maybe add QC (bcftool stats)

    emit:
    versions

}