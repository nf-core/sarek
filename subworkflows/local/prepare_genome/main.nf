//
// PREPARE GENOME
//

// Initialize channels based on params or indices that were just built
// For all modules here:
// A when clause condition is defined in the conf/modules.config to determine if the module should be run
// Condition is based on params.step and params.tools
// If and extra condition exists, it's specified in comments

include { UNTAR as UNTAR_CHR_DIR } from '../../../modules/nf-core/untar'
include { UNZIP as UNZIP_ALLELES } from '../../../modules/nf-core/unzip'
include { UNZIP as UNZIP_GC      } from '../../../modules/nf-core/unzip'
include { UNZIP as UNZIP_LOCI    } from '../../../modules/nf-core/unzip'
include { UNZIP as UNZIP_RT      } from '../../../modules/nf-core/unzip'

workflow PREPARE_GENOME {
    take:
    ascat_alleles // params.ascat_alleles
    ascat_loci    // params.ascat_loci
    ascat_loci_gc // params.ascat_loci_gc
    ascat_loci_rt // params.ascat_loci_rt
    chr_dir       // params.chr_dir

    main:
    versions = Channel.empty()

    // prepare ascat and controlfreec reference files
    if (!ascat_alleles) {
        allele_files = Channel.empty()
    }
    else if (ascat_alleles.endsWith(".zip")) {
        UNZIP_ALLELES(Channel.fromPath(file(ascat_alleles)).collect().map { it -> [[id: it[0].baseName], it] })
        allele_files = UNZIP_ALLELES.out.unzipped_archive.map { it[1] }
        versions = versions.mix(UNZIP_ALLELES.out.versions)
    }
    else {
        allele_files = Channel.fromPath(ascat_alleles).collect()
    }

    if (!ascat_loci) {
        loci_files = Channel.empty()
    }
    else if (ascat_loci.endsWith(".zip")) {
        UNZIP_LOCI(Channel.fromPath(file(ascat_loci)).collect().map { it -> [[id: it[0].baseName], it] })
        loci_files = UNZIP_LOCI.out.unzipped_archive.map { it[1] }
        versions = versions.mix(UNZIP_LOCI.out.versions)
    }
    else {
        loci_files = Channel.fromPath(ascat_loci).collect()
    }

    if (!ascat_loci_gc) {
        gc_file = Channel.value([])
    }
    else if (ascat_loci_gc.endsWith(".zip")) {
        UNZIP_GC(Channel.fromPath(file(ascat_loci_gc)).collect().map { it -> [[id: it[0].baseName], it] })
        gc_file = UNZIP_GC.out.unzipped_archive.map { it[1] }
        versions = versions.mix(UNZIP_GC.out.versions)
    }
    else {
        gc_file = Channel.fromPath(ascat_loci_gc).collect()
    }

    if (!ascat_loci_rt) {
        rt_file = Channel.value([])
    }
    else if (ascat_loci_rt.endsWith(".zip")) {
        UNZIP_RT(Channel.fromPath(file(ascat_loci_rt)).collect().map { it -> [[id: it[0].baseName], it] })
        rt_file = UNZIP_RT.out.unzipped_archive.map { it[1] }
        versions = versions.mix(UNZIP_RT.out.versions)
    }
    else {
        rt_file = Channel.fromPath(ascat_loci_rt).collect()
    }

    if (!chr_dir) {
        chr_files = Channel.value([])
    }
    else if (chr_dir.endsWith(".tar.gz")) {
        UNTAR_CHR_DIR(Channel.fromPath(file(chr_dir)).collect().map { it -> [[id: it[0].baseName], it] })
        chr_files = UNTAR_CHR_DIR.out.untar.map { it[1] }
        versions = versions.mix(UNTAR_CHR_DIR.out.versions)
    }
    else {
        chr_files = Channel.fromPath(chr_dir).collect()
    }

    emit:
    allele_files // path: allele_files
    chr_files    // path: chr_files
    gc_file      // path: gc_file
    loci_files   // path: loci_files
    rt_file      // path: rt_file
    versions     // channel: [ versions.yml ]
}
