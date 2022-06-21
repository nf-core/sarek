include { ASCAT } from '../../modules/nf-core/modules/ascat/main'

workflow RUN_ASCAT_SOMATIC {

    take:
    cram_pair                // channel: [mandatory] [meta, normal_cram, normal_crai, tumor_cram, tumor_crai]
    allele_files             // channel: [mandatory] zip
    loci_files               // channel: [mandatory] zip
    gc_file                  // channel: [optional]  txt
    rt_file                  // channel: [optional]  txt
    intervals_bed            // channel: [optional]  bed
    ref_fasta                // channel: [optional]  fasta


    main:

    ch_versions = Channel.empty()

    UNZIP_ALLELES(allele_files)
    UNZIP_LOCI(loci_files)
    //UNZIP_GC(gc_file)
    //UNZIP_RT(rt_file)

    ASCAT(cram_pair, allele_files, loci_files, [], [], [], [])

    ch_versions = ch_versions.mix(RUN_ASCAT_SOMATIC.out.versions)


    emit:
    versions = ch_versions
}
