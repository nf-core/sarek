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
include { GATK4_CREATESEQUENCEDICTIONARY         } from '../../modules/nf-core/modules/gatk4/createsequencedictionary/main'
include { MSISENSORPRO_SCAN                      } from '../../modules/nf-core/modules/msisensorpro/scan/main'
include { SAMTOOLS_FAIDX                         } from '../../modules/nf-core/modules/samtools/faidx/main'
include { TABIX_TABIX as TABIX_DBSNP             } from '../../modules/nf-core/modules/tabix/tabix/main'
include { TABIX_TABIX as TABIX_GERMLINE_RESOURCE } from '../../modules/nf-core/modules/tabix/tabix/main'
include { TABIX_TABIX as TABIX_KNOWN_INDELS      } from '../../modules/nf-core/modules/tabix/tabix/main'
include { TABIX_TABIX as TABIX_PON               } from '../../modules/nf-core/modules/tabix/tabix/main'

workflow PREPARE_GENOME {
    take:
        dbsnp             // channel: [optional]  dbsnp
        fasta             // channel: [mandatory] fasta
        fasta_fai         // channel: [optional]  fasta_fai
        germline_resource // channel: [optional]  germline_resource
        known_indels      // channel: [optional]  known_indels
        pon               // channel: [optional]  pon

    main:

    ch_versions = Channel.empty()

    BWAMEM1_INDEX(fasta) // If aligner is bwa-mem
    BWAMEM2_INDEX(fasta) // If aligner is bwa-mem2
    // if we use mix here, bwa becomes a channel that is comsumed
    ch_bwa = params.aligner == "bwa-mem" ? BWAMEM1_INDEX.out.index : BWAMEM2_INDEX.out.index

    GATK4_CREATESEQUENCEDICTIONARY(fasta)
    MSISENSORPRO_SCAN(fasta.map{ it -> [[id:it[0].baseName], it] })
    SAMTOOLS_FAIDX(fasta.map{ it -> [[id:it[0].getName()], it] })
    TABIX_DBSNP(dbsnp.map{ it -> [[id:it[0].baseName], it] })
    TABIX_GERMLINE_RESOURCE(germline_resource.map{ it -> [[id:it[0].baseName], it] })
    TABIX_KNOWN_INDELS(known_indels.map{ it -> [[id:it[0].baseName], it] })
    TABIX_PON(pon.map{ it -> [[id:it[0].baseName], it] })

    // Gather versions of all tools used
    ch_versions  = ch_versions.mix(SAMTOOLS_FAIDX.out.versions)
    ch_versions = ch_versions.mix(BWAMEM1_INDEX.out.versions)
    ch_versions = ch_versions.mix(BWAMEM2_INDEX.out.versions)
    ch_versions = ch_versions.mix(GATK4_CREATESEQUENCEDICTIONARY.out.versions)
    ch_versions = ch_versions.mix(MSISENSORPRO_SCAN.out.versions)
    ch_versions = ch_versions.mix(TABIX_DBSNP.out.versions)
    ch_versions = ch_versions.mix(TABIX_GERMLINE_RESOURCE.out.versions)
    ch_versions = ch_versions.mix(TABIX_KNOWN_INDELS.out.versions)
    ch_versions = ch_versions.mix(TABIX_PON.out.versions)

    emit:
        bwa                              = ch_bwa                                                         // path: {bwamem1,bwamem2}/index
        dbsnp_tbi                        = TABIX_DBSNP.out.tbi.map{ meta, tbi -> [tbi] }                  // path: dbsnb.vcf.gz.tbi
        dict                             = GATK4_CREATESEQUENCEDICTIONARY.out.dict                        // path: genome.fasta.dict
        fasta_fai                        = SAMTOOLS_FAIDX.out.fai.map{ meta, fai -> [fai] }               // path: genome.fasta.fai
        germline_resource_tbi            = TABIX_GERMLINE_RESOURCE.out.tbi.map{ meta, tbi -> [tbi] }      // path: germline_resource.vcf.gz.tbi
        known_indels_tbi                 = TABIX_KNOWN_INDELS.out.tbi.map{ meta, tbi -> [tbi] }.collect() // path: {known_indels*}.vcf.gz.tbi
        msisensorpro_scan                = MSISENSORPRO_SCAN.out.list.map{ meta, list -> [list] }         // path: genome_msi.list
        pon_tbi                          = TABIX_PON.out.tbi.map{ meta, tbi -> [tbi] }                    // path: pon.vcf.gz.tbi

        versions                         = ch_versions                                                    // channel: [ versions.yml ]
}
