include { BWA_INDEX as BWAMEM1_INDEX                } from '../../../modules/nf-core/bwa/index'
include { BWAMEM2_INDEX                             } from '../../../modules/nf-core/bwamem2/index'
include { DRAGMAP_HASHTABLE                         } from '../../../modules/nf-core/dragmap/hashtable'
include { GATK4_CREATESEQUENCEDICTIONARY            } from '../../../modules/nf-core/gatk4/createsequencedictionary'
include { MSISENSORPRO_SCAN                         } from '../../../modules/nf-core/msisensorpro/scan'
include { SAMTOOLS_FAIDX                            } from '../../../modules/nf-core/samtools/faidx'
include { TABIX_TABIX as TABIX_BCFTOOLS_ANNOTATIONS } from '../../../modules/nf-core/tabix/tabix'
include { TABIX_TABIX as TABIX_DBSNP                } from '../../../modules/nf-core/tabix/tabix'
include { TABIX_TABIX as TABIX_GERMLINE_RESOURCE    } from '../../../modules/nf-core/tabix/tabix'
include { TABIX_TABIX as TABIX_KNOWN_INDELS         } from '../../../modules/nf-core/tabix/tabix'
include { TABIX_TABIX as TABIX_KNOWN_SNPS           } from '../../../modules/nf-core/tabix/tabix'
include { TABIX_TABIX as TABIX_PON                  } from '../../../modules/nf-core/tabix/tabix'
include { UNTAR as UNTAR_CHR_DIR                    } from '../../../modules/nf-core/untar'
include { UNZIP as UNZIP_ALLELES                    } from '../../../modules/nf-core/unzip'
include { UNZIP as UNZIP_GC                         } from '../../../modules/nf-core/unzip'
include { UNZIP as UNZIP_LOCI                       } from '../../../modules/nf-core/unzip'
include { UNZIP as UNZIP_RT                         } from '../../../modules/nf-core/unzip'

workflow PREPARE_GENOME {
    take:
    ascat_alleles_in            // params.ascat_alleles
    ascat_loci_in               // params.ascat_loci
    ascat_loci_gc_in            // params.ascat_loci_gc
    ascat_loci_rt_in            // params.ascat_loci_rt
    bcftools_annotations_in     // params.bcftools_annotations
    bcftools_annotations_tbi_in // params.bcftools_annotations
    bwa_in                      // params.bwa
    bwamem2_in                  // params.bwamem2
    chr_dir_in                  // params.chr_dir
    dbsnp_in                    // params.dbsnp
    dbsnp_tbi_in                // params.dbsnp_tbi
    dict_in                     // params.dict
    dragmap_in                  // params.dragmap
    fasta_in                    // params.fasta
    fasta_fai_in                // params.fasta_fai
    germline_resource_in        // params.germline_resource
    germline_resource_tbi_in    // params.germline_resource_tbi
    known_indels_in             // params.known_indels
    known_indels_tbi_in         // params.known_indels_tbi
    known_snps_in               // params.known_snps
    known_snps_tbi_in           // params.known_snps_tbi
    pon_in                      // params.pon
    pon_tbi_in                  // params.pon_tbi
    aligner                     // params.aligner
    tools                       // params.tools
    step                        // params.step
    vep_include_fasta           // params.vep_include_fasta

    main:
    versions = Channel.empty()

    // TODO: EXTRACT FASTA FILE?
    fasta = fasta_in ? Channel.fromPath(fasta_in).map { fasta -> [[id: fasta.baseName], fasta] }.collect() : Channel.empty()
    vep_fasta = vep_include_fasta ? fasta : [[id: 'null'], []]

    if (!step == 'mapping') {
        index_alignment = Channel.empty()
    }
    else if (!bwa_in && (aligner == "bwa-mem" || aligner == "sentieon-bwamem" || aligner == "parabricks")) {
        BWAMEM1_INDEX(fasta)
        index_alignment = BWAMEM1_INDEX.out.index.collect()
        versions = versions.mix(BWAMEM1_INDEX.out.versions)
    }
    else if (bwa_in && (aligner == "bwa-mem" || aligner == "sentieon-bwamem" || aligner == "parabricks")) {
        index_alignment = Channel.fromPath(bwa_in).map { index -> [[id: 'bwa'], index] }.collect()
    }
    else if (!bwamem2_in && aligner == 'bwa-mem2') {
        BWAMEM2_INDEX(fasta)
        index_alignment = BWAMEM2_INDEX.out.index.collect()
        versions = versions.mix(BWAMEM2_INDEX.out.versions)
    }
    else if (bwamem2_in && aligner == 'bwa-mem2') {
        index_alignment = Channel.fromPath(bwamem2_in).map { index -> [[id: 'bwamem2'], index] }.collect()
    }
    else if (!dragmap_in && aligner == 'dragmap') {
        DRAGMAP_HASHTABLE(fasta)
        index_alignment = DRAGMAP_HASHTABLE.out.hashmap.collect()
        versions = versions.mix(DRAGMAP_HASHTABLE.out.versions)
    }
    else if (dragmap_in && aligner == 'dragmap') {
        index_alignment = Channel.fromPath(dragmap_in).map { index -> [[id: 'dragmap'], index] }.collect()
    }

    if (!dict_in && step != "annotate") {
        GATK4_CREATESEQUENCEDICTIONARY(fasta)
        dict = GATK4_CREATESEQUENCEDICTIONARY.out.dict.collect()
        versions = versions.mix(GATK4_CREATESEQUENCEDICTIONARY.out.versions)
    }
    else if (dict_in) {
        dict = Channel.fromPath(dict_in).map { it -> [[id: 'dict'], it] }.collect()
    }
    else {
        dict = Channel.empty()
    }

    if (!fasta_fai_in && step != "annotate") {
        SAMTOOLS_FAIDX(fasta, [[id: 'no_fai'], []], false)
        fasta_fai = SAMTOOLS_FAIDX.out.fai.collect()
        versions = versions.mix(SAMTOOLS_FAIDX.out.versions)
    }
    else if (fasta_fai_in) {
        fasta_fai = Channel.fromPath(fasta_fai_in).map { it -> [[id: 'fai'], it] }.collect()
    }
    else {
        fasta_fai = Channel.empty()
    }

    MSISENSORPRO_SCAN(fasta)

    if (!bcftools_annotations_in) {
        bcftools_annotations = Channel.empty()
        bcftools_annotations_tbi = Channel.empty()
    }
    else {
        bcftools_annotations = Channel.fromPath(bcftools_annotations_in).collect()
    }

    if (!bcftools_annotations_tbi_in && bcftools_annotations_in) {
        TABIX_BCFTOOLS_ANNOTATIONS(bcftools_annotations)
        bcftools_annotations_tbi_in = TABIX_BCFTOOLS_ANNOTATIONS.out.tbi.map { _meta, tbi -> [tbi] }.collect()
        versions = versions.mix(TABIX_BCFTOOLS_ANNOTATIONS.out.versions)
    }
    else {
        bcftools_annotations_tbi = Channel.fromPath(bcftools_annotations_tbi_in).collect()
    }

    if (!dbsnp_in) {
        dbsnp = Channel.empty()
        dbsnp_tbi = Channel.empty()
    }
    else {
        dbsnp = Channel.fromPath(dbsnp_in).collect()
    }

    if (!dbsnp_tbi_in && dbsnp_in) {
        TABIX_DBSNP(dbsnp)
        dbsnp_tbi_in = TABIX_DBSNP.out.tbi.map { _meta, tbi -> [tbi] }.collect()
        versions = versions.mix(TABIX_DBSNP.out.versions)
    }
    else if (dbsnp_in) {
        dbsnp_tbi = Channel.fromPath(dbsnp_tbi_in).collect()
    }
    else {
        dbsnp_tbi = Channel.empty()
    }

    if (!germline_resource_in) {
        germline_resource = Channel.empty()
        germline_resource_tbi = Channel.empty()
    }
    else {
        germline_resource = Channel.fromPath(germline_resource_in).collect()
    }

    if (!germline_resource_tbi_in && germline_resource_in) {
        TABIX_GERMLINE_RESOURCE(germline_resource)
        germline_resource_tbi_in = TABIX_GERMLINE_RESOURCE.out.tbi.map { _meta, tbi -> [tbi] }.collect()
        versions = versions.mix(TABIX_GERMLINE_RESOURCE.out.versions)
    }
    else if (germline_resource_in) {
        germline_resource_tbi = Channel.fromPath(germline_resource_tbi_in).collect()
    }
    else {
        germline_resource_tbi = Channel.empty()
    }

    if (!known_indels_in) {
        known_indels = Channel.empty()
        known_indels_tbi = Channel.empty()
    }
    else {
        known_indels = Channel.fromPath(known_indels_in).collect()
    }

    if (!known_indels_tbi_in && known_indels_in) {
        TABIX_KNOWN_INDELS(known_indels)
        known_indels_tbi_in = TABIX_KNOWN_INDELS.out.tbi.map { _meta, tbi -> [tbi] }.collect()
        versions = versions.mix(TABIX_KNOWN_INDELS.out.versions)
    }
    else if (known_indels_in) {
        known_indels_tbi = Channel.fromPath(known_indels_tbi_in).collect()
    }
    else {
        known_indels_tbi = Channel.empty()
    }

    if (!known_snps_in) {
        known_snps = Channel.empty()
        known_snps_tbi = Channel.empty()
    }
    else {
        known_snps = Channel.fromPath(known_snps_in).collect()
    }

    if (!known_snps_tbi_in && known_snps_in) {
        TABIX_KNOWN_SNPS(known_snps)
        known_snps_tbi_in = TABIX_KNOWN_SNPS.out.tbi.map { _meta, tbi -> [tbi] }.collect()
        versions = versions.mix(TABIX_KNOWN_SNPS.out.versions)
    }
    else if (known_snps_in) {
        known_snps_tbi = Channel.fromPath(known_snps_tbi_in).collect()
    }
    else {
        known_snps_tbi = Channel.empty()
    }

    if (!pon_in) {
        pon = Channel.empty()
        pon_tbi = Channel.empty()
    }
    else {
        pon = Channel.fromPath(pon_in).collect()
    }

    if (!pon_tbi_in && pon_in) {
        TABIX_PON(pon)
        pon_tbi_in = TABIX_PON.out.tbi.map { _meta, tbi -> [tbi] }.collect()
        versions = versions.mix(TABIX_PON.out.versions)
    }
    else if (pon_in) {
        pon_tbi = Channel.fromPath(pon_tbi_in).collect()
    }
    else {
        pon_tbi = Channel.empty()
    }

    // known_sites is made by grouping both the dbsnp and the known snps/indels resources
    // Which can either or both be optional
    known_sites_indels = dbsnp.concat(known_indels).collect()
    known_sites_indels_tbi = dbsnp_tbi.concat(known_indels_tbi).collect()
    known_sites_snps = dbsnp.concat(known_snps).collect()
    known_sites_snps_tbi = dbsnp_tbi.concat(known_snps_tbi).collect()

    // prepare ascat and controlfreec reference files
    if (!ascat_alleles_in) {
        ascat_alleles = Channel.empty()
    }
    else if (ascat_alleles_in.endsWith(".zip") && tools.split(',').contains('ascat')) {
        UNZIP_ALLELES(Channel.fromPath(file(ascat_alleles_in)).collect().map { it -> [[id: it[0].baseName], it] })

        ascat_alleles = UNZIP_ALLELES.out.unzipped_archive.map { _meta, extracted_archive -> extracted_archive }
        versions = versions.mix(UNZIP_ALLELES.out.versions)
    }
    else {
        ascat_alleles = Channel.fromPath(ascat_alleles_in).collect()
    }

    if (!ascat_loci_in) {
        ascat_loci = Channel.empty()
    }
    else if (ascat_loci_in.endsWith(".zip") && tools.split(',').contains('ascat')) {
        UNZIP_LOCI(Channel.fromPath(file(ascat_loci_in)).collect().map { it -> [[id: it[0].baseName], it] })

        ascat_loci = UNZIP_LOCI.out.unzipped_archive.map { _meta, extracted_archive -> extracted_archive }
        versions = versions.mix(UNZIP_LOCI.out.versions)
    }
    else {
        ascat_loci = Channel.fromPath(ascat_loci_in).collect()
    }

    if (!ascat_loci_gc_in) {
        ascat_loci_gc = Channel.value([])
    }
    else if (ascat_loci_gc_in.endsWith(".zip") && tools.split(',').contains('ascat')) {
        UNZIP_GC(Channel.fromPath(file(ascat_loci_gc_in)).collect().map { it -> [[id: it[0].baseName], it] })

        ascat_loci_gc = UNZIP_GC.out.unzipped_archive.map { _meta, extracted_archive -> extracted_archive }
        versions = versions.mix(UNZIP_GC.out.versions)
    }
    else {
        ascat_loci_gc = Channel.fromPath(ascat_loci_gc_in).collect()
    }

    if (!ascat_loci_rt_in) {
        ascat_loci_rt = Channel.value([])
    }
    else if (ascat_loci_rt_in.endsWith(".zip") && tools.split(',').contains('ascat')) {
        UNZIP_RT(Channel.fromPath(file(ascat_loci_rt_in)).collect().map { it -> [[id: it[0].baseName], it] })

        ascat_loci_rt = UNZIP_RT.out.unzipped_archive.map { _meta, extracted_archive -> extracted_archive }
        versions = versions.mix(UNZIP_RT.out.versions)
    }
    else {
        ascat_loci_rt = Channel.fromPath(ascat_loci_rt_in).collect()
    }

    if (!chr_dir_in) {
        chr_dir = Channel.value([])
    }
    else if (chr_dir_in.endsWith(".tar.gz") && tools.split(',').contains('controlfreec')) {
        UNTAR_CHR_DIR(Channel.fromPath(file(chr_dir_in)).collect().map { it -> [[id: it[0].baseName], it] })

        chr_dir = UNTAR_CHR_DIR.out.untar.map { _meta, extracted_archive -> extracted_archive }
        versions = versions.mix(UNTAR_CHR_DIR.out.versions)
    }
    else {
        chr_dir = Channel.fromPath(chr_dir_in).collect()
    }

    // Gather versions of all tools used
    versions = versions.mix(MSISENSORPRO_SCAN.out.versions)

    emit:
    ascat_alleles            // Channel: [ascat_alleles]
    ascat_loci               // Channel: [ascat_loci]
    ascat_loci_gc            // Channel: [ascat_loci_gc]
    ascat_loci_rt            // Channel: [ascat_loci_rt]
    bcftools_annotations     // Channel: [bcftools_annotations]
    bcftools_annotations_tbi // Channel: [bcftools_annotations_tbi]
    chr_dir                  // Channel: [chr_dir]
    dbsnp                    // Channel: [dbsnp]
    dbsnp_tbi                // Channel: [dbsnp_tbi]
    dict                     // Channel: [meta, dict]
    fasta                    // Channel: [meta, fasta]
    fasta_fai                // Channel: [meta, fasta_fai]
    germline_resource        // Channel: [germline_resource]
    germline_resource_tbi    // Channel: [germline_resource_tbi]
    index_alignment          // Channel: [meta, index_alignment] either bwa, bwamem2 or dragmap
    known_indels             // Channel: [known_indels]
    known_indels_tbi         // Channel: [known_indels_tbi]
    known_sites_indels       // Channel: [known_sites_indels]
    known_sites_indels_tbi   // Channel: [known_sites_indels_tbi]
    known_sites_snps         // Channel: [known_sites_snps]
    known_sites_snps_tbi     // Channel: [known_sites_snps_tbi]
    known_snps               // Channel: [known_snps]
    known_snps_tbi           // Channel: [known_snps_tbi]
    msisensorpro_scan        = MSISENSORPRO_SCAN.out.list.map { _meta, list -> [list] } // path: genome_msi.list
    pon                      // Channel: [pon]
    pon_tbi                  // Channel: [pon_tbi]
    vep_fasta                // Channel: [meta, vep_fasta]
    versions                 // channel: [versions.yml]
}
