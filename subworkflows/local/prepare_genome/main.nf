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
    ascat_alleles        // params.ascat_alleles
    ascat_loci           // params.ascat_loci
    ascat_loci_gc        // params.ascat_loci_gc
    ascat_loci_rt        // params.ascat_loci_rt
    bcftools_annotations // channel: [optional] bcftools annotations file
    bwa_in               // params.bwa
    bwamem2_in           // params.bwamem2
    chr_dir              // params.chr_dir
    dict_in              // params.dict
    dbsnp                // channel: [optional]  dbsnp
    dragmap_in           // params.dragmap
    fasta_in             // params.fasta
    fasta_fai_in         // params.fasta_fai
    germline_resource    // channel: [optional]  germline_resource
    known_indels         // channel: [optional]  known_indels
    known_snps           // channel: [optional]  known_snps
    pon                  // channel: [optional]  pon
    aligner              // params.aligner
    step                 // params.step
    vep_include_fasta    // params.vep_include_fasta

    main:
    versions = Channel.empty()

    // TODO: EXTRACT FASTA FILE?
    fasta = fasta_in ? Channel.fromPath(fasta_in).map { fasta -> [[id: fasta.baseName], fasta] }.collect() : Channel.empty()
    vep_fasta = vep_include_fasta ? fasta : [[id: 'null'], []]

    index_alignment = Channel.empty()
    dict = Channel.empty()
    fasta_fai = Channel.empty()

    if (!bwa_in && step == 'mapping' && (aligner == "bwa-mem" || aligner == "sentieon-bwamem" || aligner == "parabricks")) {
        BWAMEM1_INDEX(fasta)
        index_alignment = BWAMEM1_INDEX.out.index.collect()
        versions = versions.mix(BWAMEM1_INDEX.out.versions)
    }
    else if (bwa_in && step == 'mapping' && (aligner == "bwa-mem" || aligner == "sentieon-bwamem" || aligner == "parabricks")) {
        index_alignment = Channel.fromPath(bwa_in).map { index -> [[id: 'bwa'], index] }.collect()
    }
    else if (!bwamem2_in && step == 'mapping' && aligner == 'bwa-mem2') {
        BWAMEM2_INDEX(fasta)
        index_alignment = BWAMEM2_INDEX.out.index.collect()
        versions = versions.mix(BWAMEM2_INDEX.out.versions)
    }
    else if (bwamem2_in && step == 'mapping' && aligner == 'bwa-mem2') {
        index_alignment = Channel.fromPath(bwamem2_in).map { index -> [[id: 'bwamem2'], index] }.collect()
    }
    else if (!dragmap_in && step == 'mapping' && aligner == 'dragmap') {
        DRAGMAP_HASHTABLE(fasta)
        index_alignment = DRAGMAP_HASHTABLE.out.hashmap.collect()
        versions = versions.mix(DRAGMAP_HASHTABLE.out.versions)
    }
    else if (dragmap_in && step == 'mapping' && aligner == 'dragmap') {
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

    if (!fasta_fai_in && step != "annotate") {
        SAMTOOLS_FAIDX(fasta, [[id: 'no_fai'], []], false)
        fasta_fai = SAMTOOLS_FAIDX.out.fai.collect()
        versions = versions.mix(SAMTOOLS_FAIDX.out.versions)
    }
    else if (fasta_fai_in) {
        fasta_fai = Channel.fromPath(fasta_fai_in).map { it -> [[id: 'fai'], it] }.collect()
    }

    MSISENSORPRO_SCAN(fasta)

    // the following are flattened and mapped in case the user supplies more than one value for the param
    // written for KNOWN_INDELS, but preemptively applied to the rest
    // [ file1, file2 ] becomes [ [ meta1, file1 ], [ meta2, file2 ] ]
    // outputs are collected to maintain a single channel for relevant TBI files
    TABIX_BCFTOOLS_ANNOTATIONS(bcftools_annotations.flatten().map { it -> [[id: it.baseName], it] })
    TABIX_DBSNP(dbsnp.flatten().map { it -> [[id: it.baseName], it] })
    TABIX_GERMLINE_RESOURCE(germline_resource.flatten().map { it -> [[id: it.baseName], it] })
    TABIX_KNOWN_SNPS(known_snps.flatten().map { it -> [[id: it.baseName], it] })
    TABIX_KNOWN_INDELS(known_indels.flatten().map { it -> [[id: it.baseName], it] })
    TABIX_PON(pon.flatten().map { it -> [[id: it.baseName], it] })

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

    // Gather versions of all tools used
    versions = versions.mix(MSISENSORPRO_SCAN.out.versions)
    versions = versions.mix(TABIX_BCFTOOLS_ANNOTATIONS.out.versions)
    versions = versions.mix(TABIX_DBSNP.out.versions)
    versions = versions.mix(TABIX_GERMLINE_RESOURCE.out.versions)
    versions = versions.mix(TABIX_KNOWN_INDELS.out.versions)
    versions = versions.mix(TABIX_KNOWN_SNPS.out.versions)
    versions = versions.mix(TABIX_PON.out.versions)

    emit:
    fasta                    // Channel: [ meta, fasta ]
    index_alignment          // Channel: [ meta, index_alignment ] either bwa, bwamem2 or dragmap
    vep_fasta                // Channel: [ meta, vep_fasta ]
    dict                     // Channel: [ meta, dict ]
    fasta_fai                // Channel: [ meta, fasta_fai ]
    bcftools_annotations_tbi = TABIX_BCFTOOLS_ANNOTATIONS.out.tbi.map { meta, tbi -> [tbi] }.collect() // path: bcftools_annotations.vcf.gz.tbi
    dbsnp_tbi                = TABIX_DBSNP.out.tbi.map { meta, tbi -> [tbi] }.collect() // path: dbsnb.vcf.gz.tbi
    germline_resource_tbi    = TABIX_GERMLINE_RESOURCE.out.tbi.map { meta, tbi -> [tbi] }.collect() // path: germline_resource.vcf.gz.tbi
    known_snps_tbi           = TABIX_KNOWN_SNPS.out.tbi.map { meta, tbi -> [tbi] }.collect() // path: {known_indels*}.vcf.gz.tbi
    known_indels_tbi         = TABIX_KNOWN_INDELS.out.tbi.map { meta, tbi -> [tbi] }.collect() // path: {known_indels*}.vcf.gz.tbi
    msisensorpro_scan        = MSISENSORPRO_SCAN.out.list.map { meta, list -> [list] } // path: genome_msi.list
    pon_tbi                  = TABIX_PON.out.tbi.map { meta, tbi -> [tbi] }.collect() // path: pon.vcf.gz.tbi
    allele_files             // path: allele_files
    chr_files                // path: chr_files
    gc_file                  // path: gc_file
    loci_files               // path: loci_files
    rt_file                  // path: rt_file
    versions                 // channel: [ versions.yml ]
}
