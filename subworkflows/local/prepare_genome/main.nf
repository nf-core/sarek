//
// PREPARE GENOME
//

// Initialize channels based on params or indices that were just built
// For all modules here:
// A when clause condition is defined in the conf/modules.config to determine if the module should be run
// Condition is based on params.step and params.tools
// If and extra condition exists, it's specified in comments

include { BWA_INDEX as BWAMEM1_INDEX             } from '../../../modules/nf-core/bwa/index/main'
include { BWAMEM2_INDEX                          } from '../../../modules/nf-core/bwamem2/index/main'
include { DRAGMAP_HASHTABLE                      } from '../../../modules/nf-core/dragmap/hashtable/main'
include { GATK4_CREATESEQUENCEDICTIONARY         } from '../../../modules/nf-core/gatk4/createsequencedictionary/main'
include { MSISENSORPRO_SCAN                      } from '../../../modules/nf-core/msisensorpro/scan/main'
include { SAMTOOLS_FAIDX                         } from '../../../modules/nf-core/samtools/faidx/main'
include { TABIX_TABIX as TABIX_DBSNP             } from '../../../modules/nf-core/tabix/tabix/main'
include { TABIX_TABIX as TABIX_GERMLINE_RESOURCE } from '../../../modules/nf-core/tabix/tabix/main'
include { TABIX_TABIX as TABIX_KNOWN_INDELS      } from '../../../modules/nf-core/tabix/tabix/main'
include { TABIX_TABIX as TABIX_KNOWN_SNPS        } from '../../../modules/nf-core/tabix/tabix/main'
include { TABIX_TABIX as TABIX_PON               } from '../../../modules/nf-core/tabix/tabix/main'
include { UNTAR as UNTAR_CHR_DIR                 } from '../../../modules/nf-core/untar/main'
include { UNZIP as UNZIP_ALLELES                 } from '../../../modules/nf-core/unzip/main'
include { UNZIP as UNZIP_LOCI                    } from '../../../modules/nf-core/unzip/main'
include { UNZIP as UNZIP_GC                      } from '../../../modules/nf-core/unzip/main'
include { UNZIP as UNZIP_RT                      } from '../../../modules/nf-core/unzip/main'

workflow PREPARE_GENOME {
    take:
        ascat_alleles           // channel: [optional]  ascat allele files
        ascat_loci              // channel: [optional]  ascat loci files
        ascat_loci_gc           // channel: [optional]  ascat gc content file
        ascat_loci_rt           // channel: [optional]  ascat replictiming file
        chr_dir                 // channel: [optional]  chromosome files
        dbsnp                   // channel: [optional]  dbsnp
        fasta                   // channel: [mandatory] fasta
        fasta_fai               // channel: [optional]  fasta_fai
        germline_resource       // channel: [optional]  germline_resource
        known_indels            // channel: [optional]  known_indels
        known_snps
        pon                     // channel: [optional]  pon


    main:

    ch_versions = Channel.empty()

    BWAMEM1_INDEX(fasta.map{ it -> [[id:it[0].baseName], it] }) // If aligner is bwa-mem
    BWAMEM2_INDEX(fasta.map{ it -> [[id:it[0].baseName], it] }) // If aligner is bwa-mem2
    DRAGMAP_HASHTABLE(fasta.map{ it -> [[id:it[0].baseName], it] }) // If aligner is dragmap

    GATK4_CREATESEQUENCEDICTIONARY(fasta)
    MSISENSORPRO_SCAN(fasta.map{ it -> [[id:it[0].baseName], it] })
    SAMTOOLS_FAIDX(fasta.map{ it -> [[id:it[0].baseName], it] })

    // the following are flattened and mapped in case the user supplies more than one value for the param
    // written for KNOWN_INDELS, but preemptively applied to the rest
    // [file1, file2] becomes [[meta1, file1],[meta2, file2]]
    // outputs are collected to maintain a single channel for relevant TBI files
    TABIX_DBSNP(dbsnp.flatten().map{ it -> [[id:it.baseName], it] })
    TABIX_GERMLINE_RESOURCE(germline_resource.flatten().map{ it -> [[id:it.baseName], it] })
    TABIX_KNOWN_SNPS( known_snps.flatten().map{ it -> [[id:it.baseName], it] } )
    TABIX_KNOWN_INDELS( known_indels.flatten().map{ it -> [[id:it.baseName], it] } )
    TABIX_PON(pon.flatten().map{ it -> [[id:it.baseName], it] })

    // prepare ascat reference files
    allele_files = ascat_alleles
    if (params.ascat_alleles && params.ascat_alleles.endsWith('.zip')) {
        UNZIP_ALLELES(ascat_alleles.map{ it -> [[id:it[0].baseName], it] })
        allele_files = UNZIP_ALLELES.out.unzipped_archive.map{ it[1] }
        ch_versions = ch_versions.mix(UNZIP_ALLELES.out.versions)
    }

    loci_files = ascat_loci
    if (params.ascat_loci && params.ascat_loci.endsWith('.zip')) {
        UNZIP_LOCI(ascat_loci.map{ it -> [[id:it[0].baseName], it] })
        loci_files = UNZIP_LOCI.out.unzipped_archive.map{ it[1] }
        ch_versions = ch_versions.mix(UNZIP_LOCI.out.versions)
    }
    gc_file = ascat_loci_gc
    if (params.ascat_loci_gc && params.ascat_loci_gc.endsWith('.zip')) {
        UNZIP_GC(ascat_loci_gc.map{ it -> [[id:it[0].baseName], it] })
        gc_file = UNZIP_GC.out.unzipped_archive.map{ it[1] }
        ch_versions = ch_versions.mix(UNZIP_GC.out.versions)
    }
    rt_file = ascat_loci_rt
    if (params.ascat_loci_rt && params.ascat_loci_rt.endsWith('.zip')) {
        UNZIP_RT(ascat_loci_rt.map{ it -> [[id:it[0].baseName], it] })
        rt_file = UNZIP_RT.out.unzipped_archive.map{ it[1] }
        ch_versions = ch_versions.mix(UNZIP_RT.out.versions)
    }


    chr_files = chr_dir
    if (params.chr_dir && params.chr_dir.endsWith('tar.gz')) {
        UNTAR_CHR_DIR(chr_dir.map{ it -> [[id:it[0].baseName], it] })
        chr_files = UNTAR_CHR_DIR.out.untar.map{ it[1] }
        ch_versions = ch_versions.mix(UNTAR_CHR_DIR.out.versions)
    }

    // Gather versions of all tools used
    ch_versions = ch_versions.mix(SAMTOOLS_FAIDX.out.versions)
    ch_versions = ch_versions.mix(BWAMEM1_INDEX.out.versions)
    ch_versions = ch_versions.mix(BWAMEM2_INDEX.out.versions)
    ch_versions = ch_versions.mix(DRAGMAP_HASHTABLE.out.versions)
    ch_versions = ch_versions.mix(GATK4_CREATESEQUENCEDICTIONARY.out.versions)
    ch_versions = ch_versions.mix(MSISENSORPRO_SCAN.out.versions)
    ch_versions = ch_versions.mix(TABIX_DBSNP.out.versions)
    ch_versions = ch_versions.mix(TABIX_GERMLINE_RESOURCE.out.versions)
    ch_versions = ch_versions.mix(TABIX_KNOWN_SNPS.out.versions)
    ch_versions = ch_versions.mix(TABIX_KNOWN_INDELS.out.versions)
    ch_versions = ch_versions.mix(TABIX_PON.out.versions)

    emit:
        bwa                              = BWAMEM1_INDEX.out.index.map{ meta, index -> [index] }.collect()       // path: bwa/*
        bwamem2                          = BWAMEM2_INDEX.out.index.map{ meta, index -> [index] }.collect()       // path: bwamem2/*
        hashtable                        = DRAGMAP_HASHTABLE.out.hashmap.map{ meta, index -> [index] }.collect() // path: dragmap/*
        dbsnp_tbi                        = TABIX_DBSNP.out.tbi.map{ meta, tbi -> [tbi] }.collect()               // path: dbsnb.vcf.gz.tbi
        dict                             = GATK4_CREATESEQUENCEDICTIONARY.out.dict                               // path: genome.fasta.dict
        fasta_fai                        = SAMTOOLS_FAIDX.out.fai.map{ meta, fai -> [fai] }                      // path: genome.fasta.fai
        germline_resource_tbi            = TABIX_GERMLINE_RESOURCE.out.tbi.map{ meta, tbi -> [tbi] }.collect()   // path: germline_resource.vcf.gz.tbi
        known_snps_tbi                   = TABIX_KNOWN_SNPS.out.tbi.map{ meta, tbi -> [tbi] }.collect()          // path: {known_indels*}.vcf.gz.tbi
        known_indels_tbi                 = TABIX_KNOWN_INDELS.out.tbi.map{ meta, tbi -> [tbi] }.collect()        // path: {known_indels*}.vcf.gz.tbi
        msisensorpro_scan                = MSISENSORPRO_SCAN.out.list.map{ meta, list -> [list] }                // path: genome_msi.list
        pon_tbi                          = TABIX_PON.out.tbi.map{ meta, tbi -> [tbi] }.collect()                 // path: pon.vcf.gz.tbi
        chr_files                        = chr_files
        allele_files                     = allele_files
        loci_files                       = loci_files
        gc_file                          = gc_file
        rt_file                          = rt_file
        versions                         = ch_versions                                                         // channel: [ versions.yml ]
}
