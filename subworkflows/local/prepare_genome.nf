//
// PREPARE GENOME
//

// Initialize channels based on params or indices that were just built
// For all modules here:
// A when clause condition is defined in the conf/modules.config to determine if the module should be run
// Condition is based on params.step and params.tools
// If and extra condition exists, it's specified in comments

include { BWA_INDEX as BWAMEM1_INDEX             } from '../../modules/nf-core/modules/bwa/index/main'
include { BWAMEM2_INDEX                          } from '../../modules/nf-core/modules/bwamem2/index/main'
include { DRAGMAP_HASHTABLE                      } from '../../modules/nf-core/modules/dragmap/hashtable/main'
include { GATK4_CREATESEQUENCEDICTIONARY         } from '../../modules/nf-core/modules/gatk4/createsequencedictionary/main'
include { MSISENSORPRO_SCAN                      } from '../../modules/nf-core/modules/msisensorpro/scan/main'
include { SAMTOOLS_FAIDX                         } from '../../modules/nf-core/modules/samtools/faidx/main'
include { TABIX_TABIX as TABIX_DBSNP             } from '../../modules/nf-core/modules/tabix/tabix/main'
include { TABIX_TABIX as TABIX_GERMLINE_RESOURCE } from '../../modules/nf-core/modules/tabix/tabix/main'
include { TABIX_TABIX as TABIX_KNOWN_INDELS      } from '../../modules/nf-core/modules/tabix/tabix/main'
include { TABIX_TABIX as TABIX_PON               } from '../../modules/nf-core/modules/tabix/tabix/main'
include { UNTAR as UNTAR_CHR_DIR                 } from '../../modules/nf-core/modules/untar/main'

workflow PREPARE_GENOME {
    take:
        chr_dir           // channel: [optional]  chromosome files
        dbsnp             // channel: [optional]  dbsnp
        fasta             // channel: [mandatory] fasta
        fasta_fai         // channel: [optional]  fasta_fai
        germline_resource // channel: [optional]  germline_resource
        known_indels      // channel: [optional]  known_indels
        pon               // channel: [optional]  pon

    main:

    ch_versions = Channel.empty()

    BWAMEM1_INDEX(fasta)     // If aligner is bwa-mem
    BWAMEM2_INDEX(fasta)     // If aligner is bwa-mem2
    DRAGMAP_HASHTABLE(fasta) // If aligner is dragmap

    GATK4_CREATESEQUENCEDICTIONARY(fasta)
    MSISENSORPRO_SCAN(fasta.map{ it -> [[id:it[0].baseName], it] })
    SAMTOOLS_FAIDX(fasta.map{ it -> [[id:it[0].getName()], it] })

    // the following are flattened and mapped in case the user supplies more than one value for the param
    // written for KNOWN_INDELS, but preemptively applied to the rest
    // [file1,file2] becomes [[meta1,file1],[meta2,file2]]
    // outputs are collected to maintain a single channel for relevant TBI files
    TABIX_DBSNP(dbsnp.flatten().map{ it -> [[id:it.baseName], it] })
    TABIX_GERMLINE_RESOURCE(germline_resource.flatten().map{ it -> [[id:it.baseName], it] })
    TABIX_KNOWN_INDELS( known_indels.flatten().map{ it -> [[id:it.baseName], it] } )
    TABIX_PON(pon.flatten().map{ it -> [[id:it.baseName], it] })

    chr_files = chr_dir
    if ( params.chr_dir.endsWith('tar.gz')){
        UNTAR_CHR_DIR(chr_dir.map{ it -> [[id:it[0].baseName], it] })
        chr_files = UNTAR_CHR_DIR.out.untar.map{ it[1] }
        ch_versions = ch_versions.mix(UNTAR_CHR_DIR.out.versions)
    }

    // Gather versions of all tools used
    ch_versions = ch_versions.mix(SAMTOOLS_FAIDX.out.versions)
    ch_versions = ch_versions.mix(BWAMEM1_INDEX.out.versions)
    ch_versions = ch_versions.mix(BWAMEM2_INDEX.out.versions)
    ch_versions = ch_versions.mix(GATK4_CREATESEQUENCEDICTIONARY.out.versions)
    ch_versions = ch_versions.mix(MSISENSORPRO_SCAN.out.versions)
    ch_versions = ch_versions.mix(TABIX_DBSNP.out.versions)
    ch_versions = ch_versions.mix(TABIX_GERMLINE_RESOURCE.out.versions)
    ch_versions = ch_versions.mix(TABIX_KNOWN_INDELS.out.versions)
    ch_versions = ch_versions.mix(TABIX_PON.out.versions)

    emit:
        bwa                              = BWAMEM1_INDEX.out.index                                              // path: bwa/*
        bwamem2                          = BWAMEM2_INDEX.out.index                                              // path: bwamem2/*
        hashtable                        = DRAGMAP_HASHTABLE.out.hashmap                                        // path: dragmap/*
        dbsnp_tbi                        = TABIX_DBSNP.out.tbi.map{ meta, tbi -> [tbi] }.collect()              // path: dbsnb.vcf.gz.tbi
        dict                             = GATK4_CREATESEQUENCEDICTIONARY.out.dict                              // path: genome.fasta.dict
        fasta_fai                        = SAMTOOLS_FAIDX.out.fai.map{ meta, fai -> [fai] }                     // path: genome.fasta.fai
        germline_resource_tbi            = TABIX_GERMLINE_RESOURCE.out.tbi.map{ meta, tbi -> [tbi] }.collect()  // path: germline_resource.vcf.gz.tbi
        known_indels_tbi                 = TABIX_KNOWN_INDELS.out.tbi.map{ meta, tbi -> [tbi] }.collect()                                                        // path: {known_indels*}.vcf.gz.tbi
        msisensorpro_scan                = MSISENSORPRO_SCAN.out.list.map{ meta, list -> [list] } // path: genome_msi.list
        pon_tbi                          = TABIX_PON.out.tbi.map{ meta, tbi -> [tbi] }.collect()                // path: pon.vcf.gz.tbi
        chr_files                        = chr_files
        versions                         = ch_versions                                            // channel: [ versions.yml ]
}
