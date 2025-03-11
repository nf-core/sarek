//
// PREPARE GENOME
//

// Initialize channels based on params or indices that were just built
// For all modules here:
// A when clause condition is defined in the conf/modules.config to determine if the module should be run
// Condition is based on params.step and params.tools
// If and extra condition exists, it's specified in comments

include { BWA_INDEX as BWAMEM1_INDEX                } from '../../../modules/nf-core/bwa/index/main'
include { BWAMEM2_INDEX                             } from '../../../modules/nf-core/bwamem2/index/main'
include { DRAGMAP_HASHTABLE                         } from '../../../modules/nf-core/dragmap/hashtable/main'
include { GATK4_CREATESEQUENCEDICTIONARY            } from '../../../modules/nf-core/gatk4/createsequencedictionary/main'
include { MINIMAP2_INDEX                            } from '../../../modules/nf-core/minimap2/index/main'
include { MSISENSORPRO_SCAN                         } from '../../../modules/nf-core/msisensorpro/scan/main'
include { SAMTOOLS_FAIDX                            } from '../../../modules/nf-core/samtools/faidx/main'
include { TABIX_TABIX as TABIX_BCFTOOLS_ANNOTATIONS } from '../../../modules/nf-core/tabix/tabix/main'
include { TABIX_TABIX as TABIX_DBSNP                } from '../../../modules/nf-core/tabix/tabix/main'
include { TABIX_TABIX as TABIX_GERMLINE_RESOURCE    } from '../../../modules/nf-core/tabix/tabix/main'
include { TABIX_TABIX as TABIX_KNOWN_INDELS         } from '../../../modules/nf-core/tabix/tabix/main'
include { TABIX_TABIX as TABIX_KNOWN_SNPS           } from '../../../modules/nf-core/tabix/tabix/main'
include { TABIX_TABIX as TABIX_PON                  } from '../../../modules/nf-core/tabix/tabix/main'
include { UNTAR as UNTAR_CHR_DIR                    } from '../../../modules/nf-core/untar/main'
include { UNZIP as UNZIP_ALLELES                    } from '../../../modules/nf-core/unzip/main'
include { UNZIP as UNZIP_GC                         } from '../../../modules/nf-core/unzip/main'
include { UNZIP as UNZIP_LOCI                       } from '../../../modules/nf-core/unzip/main'
include { UNZIP as UNZIP_RT                         } from '../../../modules/nf-core/unzip/main'

workflow PREPARE_GENOME {
    take:
    ascat_alleles        // params.ascat_alleles
    ascat_loci           // params.ascat_loci
    ascat_loci_gc        // params.ascat_loci_gc
    ascat_loci_rt        // params.ascat_loci_rt
    bcftools_annotations // channel: [optional] bcftools annotations file
    chr_dir              // params.chr_dir
    dbsnp                // channel: [optional]  dbsnp
    fasta                // channel: [mandatory] fasta
    germline_resource    // channel: [optional]  germline_resource
    known_indels         // channel: [optional]  known_indels
    known_snps           // channel: [optional]  known_snps
    pon                  // channel: [optional]  pon


    main:
    versions = Channel.empty()

    BWAMEM1_INDEX(fasta)     // If aligner is bwa-mem
    BWAMEM2_INDEX(fasta)     // If aligner is bwa-mem2
    DRAGMAP_HASHTABLE(fasta) // If aligner is dragmap
    MINIMAP2_INDEX(fasta)    // If aligner is minimap2

    GATK4_CREATESEQUENCEDICTIONARY(fasta)
    MSISENSORPRO_SCAN(fasta)
    SAMTOOLS_FAIDX(fasta, [ [ id:'no_fai' ], [] ] )

    // the following are flattened and mapped in case the user supplies more than one value for the param
    // written for KNOWN_INDELS, but preemptively applied to the rest
    // [ file1, file2 ] becomes [ [ meta1, file1 ], [ meta2, file2 ] ]
    // outputs are collected to maintain a single channel for relevant TBI files
    TABIX_BCFTOOLS_ANNOTATIONS(bcftools_annotations.flatten().map{ it -> [ [ id:it.baseName ], it ] })
    TABIX_DBSNP(dbsnp.flatten().map{ it -> [ [ id:it.baseName ], it ] })
    TABIX_GERMLINE_RESOURCE(germline_resource.flatten().map{ it -> [ [ id:it.baseName ], it ] })
    TABIX_KNOWN_SNPS(known_snps.flatten().map{ it -> [ [ id:it.baseName ], it ] } )
    TABIX_KNOWN_INDELS(known_indels.flatten().map{ it -> [ [ id:it.baseName ], it ] } )
    TABIX_PON(pon.flatten().map{ it -> [ [ id:it.baseName ], it ] })

    // prepare ascat and controlfreec reference files
    if (!ascat_alleles) allele_files = Channel.empty()
    else if (ascat_alleles.endsWith(".zip")) {
        UNZIP_ALLELES(Channel.fromPath(file(ascat_alleles)).collect().map{ it -> [ [ id:it[0].baseName ], it ] })
        allele_files = UNZIP_ALLELES.out.unzipped_archive.map{ it[1] }
        versions = versions.mix(UNZIP_ALLELES.out.versions)
    } else allele_files = Channel.fromPath(ascat_alleles).collect()

    if (!ascat_loci) loci_files = Channel.empty()
    else if (ascat_loci.endsWith(".zip")) {
        UNZIP_LOCI(Channel.fromPath(file(ascat_loci)).collect().map{ it -> [ [ id:it[0].baseName ], it ] })
        loci_files = UNZIP_LOCI.out.unzipped_archive.map{ it[1] }
        versions = versions.mix(UNZIP_LOCI.out.versions)
    } else loci_files = Channel.fromPath(ascat_loci).collect()

    if (!ascat_loci_gc) gc_file = Channel.value([])
    else if (ascat_loci_gc.endsWith(".zip")) {
        UNZIP_GC(Channel.fromPath(file(ascat_loci_gc)).collect().map{ it -> [ [ id:it[0].baseName ], it ] })
        gc_file = UNZIP_GC.out.unzipped_archive.map{ it[1] }
        versions = versions.mix(UNZIP_GC.out.versions)
    } else gc_file = Channel.fromPath(ascat_loci_gc).collect()

    if (!ascat_loci_rt) rt_file = Channel.value([])
    else if (ascat_loci_rt.endsWith(".zip")) {
        UNZIP_RT(Channel.fromPath(file(ascat_loci_rt)).collect().map{ it -> [ [ id:it[0].baseName ], it ] })
        rt_file = UNZIP_RT.out.unzipped_archive.map{ it[1] }
        versions = versions.mix(UNZIP_RT.out.versions)
    } else rt_file = Channel.fromPath(ascat_loci_rt).collect()

    if (!chr_dir) chr_files = Channel.value([])
    else if (chr_dir.endsWith(".tar.gz")) {
        UNTAR_CHR_DIR(Channel.fromPath(file(chr_dir)).collect().map{ it -> [ [ id:it[0].baseName ], it ] })
        chr_files = UNTAR_CHR_DIR.out.untar.map{ it[1] }
        versions = versions.mix(UNTAR_CHR_DIR.out.versions)
    } else chr_files = Channel.fromPath(chr_dir).collect()

    // Gather versions of all tools used
    versions = versions.mix(BWAMEM1_INDEX.out.versions)
    versions = versions.mix(BWAMEM2_INDEX.out.versions)
    versions = versions.mix(DRAGMAP_HASHTABLE.out.versions)
    versions = versions.mix(GATK4_CREATESEQUENCEDICTIONARY.out.versions)
    versions = versions.mix(MSISENSORPRO_SCAN.out.versions)
    versions = versions.mix(MINIMAP2_INDEX.out.versions)
    versions = versions.mix(SAMTOOLS_FAIDX.out.versions)
    versions = versions.mix(TABIX_BCFTOOLS_ANNOTATIONS.out.versions)
    versions = versions.mix(TABIX_DBSNP.out.versions)
    versions = versions.mix(TABIX_GERMLINE_RESOURCE.out.versions)
    versions = versions.mix(TABIX_KNOWN_INDELS.out.versions)
    versions = versions.mix(TABIX_KNOWN_SNPS.out.versions)
    versions = versions.mix(TABIX_PON.out.versions)

    emit:
    bcftools_annotations_tbi = TABIX_BCFTOOLS_ANNOTATIONS.out.tbi.map{ meta, tbi -> [tbi] }.collect()   // path: bcftools_annotations.vcf.gz.tbi
    bwa                      = BWAMEM1_INDEX.out.index.collect()                                        // path: bwa/*
    bwamem2                  = BWAMEM2_INDEX.out.index.collect()                                        // path: bwamem2/*
    hashtable                = DRAGMAP_HASHTABLE.out.hashmap.collect()                                  // path: dragmap/*
    dbsnp_tbi                = TABIX_DBSNP.out.tbi.map{ meta, tbi -> [tbi] }.collect()                  // path: dbsnb.vcf.gz.tbi
    dict                     = GATK4_CREATESEQUENCEDICTIONARY.out.dict.collect()                        // path: genome.fasta.dict
    fasta_fai                = SAMTOOLS_FAIDX.out.fai.collect()                                         // path: genome.fasta.fai
    germline_resource_tbi    = TABIX_GERMLINE_RESOURCE.out.tbi.map{ meta, tbi -> [tbi] }.collect()      // path: germline_resource.vcf.gz.tbi
    known_snps_tbi           = TABIX_KNOWN_SNPS.out.tbi.map{ meta, tbi -> [tbi] }.collect()             // path: {known_indels*}.vcf.gz.tbi
    known_indels_tbi         = TABIX_KNOWN_INDELS.out.tbi.map{ meta, tbi -> [tbi] }.collect()           // path: {known_indels*}.vcf.gz.tbi
    minimap2                 = MINIMAP2_INDEX.out.index.collect()                                       // path: genome.mmi
    msisensorpro_scan        = MSISENSORPRO_SCAN.out.list.map{ meta, list -> [list] }                   // path: genome_msi.list
    pon_tbi                  = TABIX_PON.out.tbi.map{ meta, tbi -> [tbi] }.collect()                    // path: pon.vcf.gz.tbi

    allele_files    // path: allele_files
    chr_files       // path: chr_files
    gc_file         // path: gc_file
    loci_files      // path: loci_files
    rt_file         // path: rt_file

    versions        // channel: [ versions.yml ]
}
