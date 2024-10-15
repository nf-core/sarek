
//
// Run BCFtools to annotate VCF files
//

include { BCFTOOLS_ANNOTATE } from '../../../modules/nf-core/bcftools/annotate/main'

workflow VCF_ANNOTATE_BCFTOOLS {
    take:
    vcf               // channel: [ val(meta), vcf ]
    annotations       //
    annotations_index //
    header_lines      //


    main:
    ch_versions = Channel.empty()

    BCFTOOLS_ANNOTATE(vcf.map{ meta, vcf -> [ meta, vcf, [] ] }, annotations, annotations_index, header_lines)

    ch_vcf_tbi = BCFTOOLS_ANNOTATE.out.vcf.join(BCFTOOLS_ANNOTATE.out.tbi, failOnDuplicate: true, failOnMismatch: true)

    // Gather versions of all tools used
    ch_versions = ch_versions.mix(BCFTOOLS_ANNOTATE.out.versions)

    emit:
    vcf_tbi  = ch_vcf_tbi  // channel: [ val(meta), vcf.gz, vcf.gz.tbi ]
    versions = ch_versions                  //    path: versions.yml
}
