//
// PREPARE GENOME
//

// Initialize channels based on params or indices that were just built

include { BUILD_INTERVALS                                     } from '../../modules/local/build_intervals/main'
include { BWA_INDEX as BWAMEM1_INDEX                          } from '../../modules/nf-core/modules/bwa/index/main'
include { BWAMEM2_INDEX                                       } from '../../modules/nf-core/modules/bwamem2/index/main'
include { CREATE_INTERVALS_BED                                } from '../../modules/local/create_intervals_bed/main'
include { GATK4_CREATESEQUENCEDICTIONARY                      } from '../../modules/nf-core/modules/gatk4/createsequencedictionary/main'
include { MSISENSORPRO_SCAN                                   } from '../../modules/nf-core/modules/msisensorpro/scan/main'
include { SAMTOOLS_FAIDX                                      } from '../../modules/nf-core/modules/samtools/faidx/main'
include { TABIX_BGZIPTABIX as TABIX_BGZIPTABIX_INTERVAL_SPLIT } from '../../modules/nf-core/modules/tabix/bgziptabix/main'
include { TABIX_BGZIPTABIX as TABIX_BGZIPTABIX_INTERVAL_ALL   } from '../../modules/nf-core/modules/tabix/bgziptabix/main'
include { TABIX_TABIX as TABIX_DBSNP                          } from '../../modules/nf-core/modules/tabix/tabix/main'
include { TABIX_TABIX as TABIX_GERMLINE_RESOURCE              } from '../../modules/nf-core/modules/tabix/tabix/main'
include { TABIX_TABIX as TABIX_KNOWN_INDELS                   } from '../../modules/nf-core/modules/tabix/tabix/main'
include { TABIX_TABIX as TABIX_PON                            } from '../../modules/nf-core/modules/tabix/tabix/main'

workflow PREPARE_GENOME {
    take:
        dbsnp             // channel: [optional]  dbsnp
        fasta             // channel: [mandatory] fasta
        fasta_fai         // channel: [optional]  fasta_fai
        germline_resource // channel: [optional]  germline_resource
        known_indels      // channel: [optional]  known_indels
        pon               // channel: [optional]  pon
        tools             // value:   [mandatory] tools
        step              // value:   [mandatory] step

    main:

    if(!tools) tools = ""

    ch_versions = Channel.empty()

    ch_bwa = Channel.empty()
    if (!(params.bwa) && 'mapping' in step) {
        if (params.aligner == "bwa-mem") {
            BWAMEM1_INDEX(fasta)
            ch_bwa = BWAMEM1_INDEX.out.index
            ch_versions = ch_versions.mix(BWAMEM1_INDEX.out.versions)
        } else if (params.aligner == "bwa-mem2") {
            BWAMEM2_INDEX(fasta)
            ch_bwa = BWAMEM2_INDEX.out.index
            ch_versions = ch_versions.mix(BWAMEM2_INDEX.out.versions)
        }
    }

    ch_dict = Channel.empty()
    if (!(params.dict) && !('annotate' in step) && !('controlfreec' in step)) {
        GATK4_CREATESEQUENCEDICTIONARY(fasta)
        ch_dict = GATK4_CREATESEQUENCEDICTIONARY.out.dict
        ch_versions = ch_versions.mix(GATK4_CREATESEQUENCEDICTIONARY.out.versions)
    }

    ch_fasta_fai = Channel.empty()
    if (fasta_fai) ch_fasta_fai = fasta_fai
    if (!(params.fasta_fai) && !('annotate' in step)) {
        SAMTOOLS_FAIDX(fasta.map{ it -> [[id:it[0].getName()], it]})
        ch_fasta_fai = SAMTOOLS_FAIDX.out.fai.map{ meta, fai -> [fai] }
        ch_versions = ch_versions.mix(SAMTOOLS_FAIDX.out.versions)
    }

    ch_dbsnp_tbi = Channel.empty()
    if (!(params.dbsnp_tbi) && params.dbsnp && ('mapping' in step || 'prepare_recalibration' in step || tools.contains('controlfreec') || tools.contains('haplotypecaller') || tools.contains('mutect2') || tools.contains('tnscope'))) {
        TABIX_DBSNP(dbsnp.map{ it -> [[id:it[0].baseName], it] })
        ch_dbsnp_tbi = TABIX_DBSNP.out.tbi.map{ meta, tbi -> [tbi] }
        ch_versions = ch_versions.mix(TABIX_DBSNP.out.versions)
    }

    ch_germline_resource_tbi = Channel.empty()
    if (!(params.germline_resource_tbi) && params.germline_resource && tools.contains('mutect2')) {
        TABIX_GERMLINE_RESOURCE(germline_resource.map{ it -> [[id:it[0].baseName], it] })
        ch_germline_resource_tbi = TABIX_GERMLINE_RESOURCE.out.tbi.map{ meta, tbi -> [tbi] }
        ch_versions = ch_versions.mix(TABIX_GERMLINE_RESOURCE.out.versions)
    }

    ch_known_indels_tbi = Channel.empty()
    if (!(params.known_indels_tbi) && params.known_indels && ('mapping' in step || 'prepare_recalibration' in step)) {
        TABIX_KNOWN_INDELS(known_indels.map{ it -> [[id:it[0].baseName], it] })
        ch_known_indels_tbi = TABIX_KNOWN_INDELS.out.tbi.map{ meta, tbi -> [tbi] }
        ch_versions = ch_versions.mix(TABIX_KNOWN_INDELS.out.versions)
    }

    ch_pon_tbi = Channel.empty()
    if (!(params.pon_tbi) && params.pon && (tools.contains('tnscope') || tools.contains('mutect2'))) {
        TABIX_PON(pon.map{ it -> [[id:it[0].baseName], it] })
        ch_pon_tbi = TABIX_PON.out.tbi.map{ meta, tbi -> [tbi] }
        ch_versions = ch_versions.mix(TABIX_PON.out.versions)
    }

    ch_msisensorpro_scan = Channel.empty()
    if (tools.contains('msisensorpro')) {
        MSISENSORPRO_SCAN(fasta.map{it -> [[id:it[0].baseName], it]})
        ch_msisensorpro_scan = MSISENSORPRO_SCAN.out.list.map{ meta, list -> list}
        ch_versions = ch_versions.mix(MSISENSORPRO_SCAN.out.versions)
    }

    ch_intervals                        = Channel.empty()
    ch_intervals_bed_gz_tbi             = Channel.empty()
    ch_intervals_combined_bed_gz_tbi    = Channel.empty()     //Create bed.gz and bed.gz.tbi for input/or created interval file. It contains ALL regions.

    tabix_in_combined = Channel.empty()
    if (params.no_intervals) {

        file("${params.outdir}/no_intervals.bed").text = "no_intervals\n"
        ch_intervals = Channel.fromPath(file("${params.outdir}/no_intervals.bed"))
        tabix_in_combined = ch_intervals.map{it -> [[id:it.getName()], it] }

    } else if (!('annotate' in step) && !('controlfreec' in step)) {
        if (!params.intervals){

            BUILD_INTERVALS(ch_fasta_fai)
            tabix_in_combined = BUILD_INTERVALS.out.bed.map{it -> [[id:it.getName()], it] }
            ch_intervals = CREATE_INTERVALS_BED(BUILD_INTERVALS.out.bed)

        }else{

            tabix_in_combined = Channel.fromPath(file(params.intervals)).map{it -> [[id:it.baseName], it] }
            ch_intervals = CREATE_INTERVALS_BED(file(params.intervals))

        }
    }

    if (!('annotate' in step) && !('controlfreec' in step)){
        TABIX_BGZIPTABIX_INTERVAL_ALL(tabix_in_combined)
        ch_intervals_combined_bed_gz_tbi = TABIX_BGZIPTABIX_INTERVAL_ALL.out.gz_tbi.map{ meta, bed, tbi -> [bed, tbi] }
        ch_versions = ch_versions.mix(TABIX_BGZIPTABIX_INTERVAL_ALL.out.versions)

        if (!params.no_intervals) {
            ch_intervals = ch_intervals.flatten()
                .map{ intervalFile ->
                    def duration = 0.0
                    for (line in intervalFile.readLines()) {
                        final fields = line.split('\t')
                        if (fields.size() >= 5) duration += fields[4].toFloat()
                        else {
                            start = fields[1].toInteger()
                            end = fields[2].toInteger()
                            duration += (end - start) / params.nucleotides_per_second
                        }
                    }
                    [duration, intervalFile]
                }.toSortedList({ a, b -> b[0] <=> a[0] })
                .flatten().collate(2)
                .map{duration, intervalFile -> intervalFile}
        }

        // Create bed.gz and bed.gz.tbi for each interval file. They are split by region (see above)
        tabix_in = ch_intervals.map{it -> [[id:it.baseName], it] }
        TABIX_BGZIPTABIX_INTERVAL_SPLIT(tabix_in)
        ch_intervals_bed_gz_tbi = TABIX_BGZIPTABIX_INTERVAL_SPLIT.out.gz_tbi.map{ meta, bed, tbi -> [bed, tbi] }
        ch_versions = ch_versions.mix(TABIX_BGZIPTABIX_INTERVAL_SPLIT.out.versions)
    }

    emit:
        bwa                              = ch_bwa                               // path: {bwa,bwamem2}/index
        dbsnp_tbi                        = ch_dbsnp_tbi                         // path: dbsnb.vcf.gz.tbi
        dict                             = ch_dict                              // path: genome.fasta.dict
        fasta_fai                        = ch_fasta_fai                         // path: genome.fasta.fai
        germline_resource_tbi            = ch_germline_resource_tbi             // path: germline_resource.vcf.gz.tbi
        known_indels_tbi                 = ch_known_indels_tbi.collect()        // path: {known_indels*}.vcf.gz.tbi
        msisensorpro_scan                = ch_msisensorpro_scan                 // path: genome_msi.list
        pon_tbi                          = ch_pon_tbi                           // path: pon.vcf.gz.tbi
        intervals_bed                    = ch_intervals                         // path: intervals.bed                        [intervals split for parallel execution]
        intervals_bed_gz_tbi             = ch_intervals_bed_gz_tbi              // path: target.bed.gz, target.bed.gz.tbi     [intervals split for parallel execution]
        intervals_combined_bed_gz_tbi    = ch_intervals_combined_bed_gz_tbi     // path: interval.bed.gz, interval.bed.gz.tbi [all intervals in one file]

        versions                         = ch_versions                          // channel: [ versions.yml ]
}
