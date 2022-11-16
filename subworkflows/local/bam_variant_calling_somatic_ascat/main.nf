include { ASCAT } from '../../../modules/nf-core/ascat/main'

workflow BAM_VARIANT_CALLING_SOMATIC_ASCAT {

    take:
    cram_pair                // channel: [mandatory] [meta, normal_cram, normal_crai, tumor_cram, tumor_crai]
    allele_files             // channel: [mandatory] zip
    loci_files               // channel: [mandatory] zip
    intervals_bed            // channel: [optional]  bed for WES
    fasta                    // channel: [optional]  fasta needed for cram
    gc_file                  // channel: [optional]  txt for LogRCorrection
    rt_file                  // channel: [optional]  txt for LogRCorrection

    main:

    ch_versions = Channel.empty()

    if (params.wes){
        ASCAT(cram_pair, allele_files, loci_files, intervals_bed, fasta, gc_file, rt_file)
    } else if (!params.wes) {
        ASCAT(cram_pair, allele_files, loci_files, [], fasta, gc_file, rt_file)
    }

    ch_versions = ch_versions.mix(ASCAT.out.versions)

    emit:
    versions = ch_versions
}
