//
// BUILDING INDICES
//

params.bgziptabix_target_bed_options   = [:]
params.build_intervals_options         = [:]
params.bwa_index_options               = [:]
params.bwamem2_index_options           = [:]
params.create_intervals_bed_options    = [:]
params.gatk4_dict_options              = [:]
params.msisensorpro_scan_options       = [:]
params.samtools_faidx_options          = [:]
params.tabix_dbsnp_options             = [:]
params.tabix_germline_resource_options = [:]
params.tabix_known_indels_options      = [:]
params.tabix_pon_options               = [:]

// Initialize channels based on params or indices that were just built

include { BUILD_INTERVALS }                        from '../../modules/local/build_intervals/main'                           addParams(options: params.build_intervals_options)
include { BWA_INDEX as BWAMEM1_INDEX }             from '../../modules/nf-core/modules/bwa/index/main'                       addParams(options: params.bwa_index_options)
include { BWAMEM2_INDEX }                          from '../../modules/nf-core/modules/bwamem2/index/main'                   addParams(options: params.bwamem2_index_options)
include { CREATE_INTERVALS_BED }                   from '../../modules/local/create_intervals_bed/main'                      addParams(options: params.create_intervals_bed_options)
include { GATK4_CREATESEQUENCEDICTIONARY  }        from '../../modules/nf-core/modules/gatk4/createsequencedictionary/main'  addParams(options: params.gatk4_dict_options)
include { MSISENSORPRO_SCAN }                      from '../../modules/local/msisensorpro/scan/main'                         addParams(options: params.msisensorpro_scan_options)
include { SAMTOOLS_FAIDX }                         from '../../modules/nf-core/modules/samtools/faidx/main'                  addParams(options: params.samtools_faidx_options)
include { TABIX_BGZIPTABIX }                       from '../../modules/nf-core/modules/tabix/bgziptabix/main'                addParams(options: params.bgziptabix_target_bed_options)
include { TABIX_TABIX as TABIX_DBSNP }             from '../../modules/nf-core/modules/tabix/tabix/main'                     addParams(options: params.tabix_dbsnp_options)
include { TABIX_TABIX as TABIX_GERMLINE_RESOURCE } from '../../modules/nf-core/modules/tabix/tabix/main'                     addParams(options: params.tabix_germline_resource_options)
include { TABIX_TABIX as TABIX_KNOWN_INDELS }      from '../../modules/nf-core/modules/tabix/tabix/main'                     addParams(options: params.tabix_known_indels_options)
include { TABIX_TABIX as TABIX_PON }               from '../../modules/nf-core/modules/tabix/tabix/main'                     addParams(options: params.tabix_pon_options)

workflow PREPARE_GENOME {
    take:
        dbsnp             // channel: [optional]  dbsnp
        fasta             // channel: [mandatory] fasta
        fasta_fai         // channel: [optional]  fasta_fai
        germline_resource // channel: [optional]  germline_resource
        known_indels      // channel: [optional]  known_indels
        pon               // channel: [optional]  pon
        target_bed        // channel: [optionnal] target_bed
        tools
        step

    main:

    ch_versions = Channel.empty()

    result_bwa = Channel.empty()
    if (!(params.bwa) && 'mapping' in step) {
        if (params.aligner == "bwa-mem") {
            BWAMEM1_INDEX(fasta)
            result_bwa = BWAMEM1_INDEX.out.index
            ch_versions = ch_versions.mix(BWAMEM1_INDEX.out.versions)
        } else if (params.aligner == "bwa-mem2") {
            BWAMEM2_INDEX(fasta)
            result_bwa = BWAMEM2_INDEX.out.index
            ch_versions = ch_versions.mix(BWAMEM2_INDEX.out.versions)
        }
    }

    result_dict = Channel.empty()
    if (!(params.dict) && !('annotate' in step) && !('controlfreec' in step)) {
        GATK4_CREATESEQUENCEDICTIONARY(fasta)
        result_dict = GATK4_CREATESEQUENCEDICTIONARY.out.dict
        ch_versions = ch_versions.mix(GATK4_CREATESEQUENCEDICTIONARY.out.versions)
    }

    result_fai = Channel.empty()
    if (fasta_fai) result_fai = fasta_fai
    if (!(params.fasta_fai) && !('annotate' in step)) {
        SAMTOOLS_FAIDX(fasta)
        result_fai = SAMTOOLS_FAIDX.out.fai
        ch_versions = ch_versions.mix(SAMTOOLS_FAIDX.out.versions)
    }

    ch_tabix_versions = Channel.empty()

    result_target_bed = Channel.empty()
    if ((params.target_bed) && ('manta' in tools || 'strelka' in tools)) {
        target_bed_id = target_bed.map{ it -> [[id:it[0].getName()], it] }
        TABIX_BGZIPTABIX(target_bed_id)
        result_target_bed = TABIX_BGZIPTABIX.out.tbi.map{ meta, bed, tbi -> [bed, tbi] }
        ch_tabix_versions = ch_tabix_versions.mix(TABIX_BGZIPTABIX.out.versions).first()
    }

    result_dbsnp_tbi = Channel.empty()
    if (!(params.dbsnp_tbi) && params.dbsnp && ('mapping' in step || 'preparerecalibration' in step || 'controlfreec' in tools || 'haplotypecaller' in tools|| 'mutect2' in tools || 'tnscope' in tools)) {
        dbsnp_id = dbsnp.map{ it -> [[id:it[0].baseName], it] }
        TABIX_DBSNP(dbsnp_id)
        result_dbsnp_tbi = TABIX_DBSNP.out.tbi.map{ meta, tbi -> [dbsnp, tbi] }
        ch_tabix_versions = ch_tabix_versions.mix(TABIX_DBSNP.out.versions).first()
    }

    result_germline_resource_tbi = Channel.empty()
    if (!(params.germline_resource_tbi) && params.germline_resource && 'mutect2' in tools) {
        germline_resource_id = germline_resource.map{ it -> [[id:it[0].baseName], it] }
        TABIX_GERMLINE_RESOURCE(germline_resource_id)
        result_germline_resource_tbi = TABIX_GERMLINE_RESOURCE.out.tbi.map{ meta, tbi -> [germline_resource, tbi] }
        ch_tabix_versions = ch_tabix_versions.mix(TABIX_GERMLINE_RESOURCE.out.versions).first()
    }

    result_known_indels_tbi = Channel.empty()
    if (!(params.known_indels_tbi) && params.known_indels && ('mapping' in step || 'preparerecalibration' in step)) {
        known_indels_id = known_indels.map{ it -> [[id:it[0].baseName], it] }
        TABIX_KNOWN_INDELS(known_indels_id)
        result_known_indels_tbi = TABIX_KNOWN_INDELS.out.tbi.map{ meta, tbi -> [known_indels, tbi] }
        ch_tabix_versions = ch_tabix_versions.mix(TABIX_KNOWN_INDELS.out.versions).first()
    }

    result_pon_tbi = Channel.empty()
    if (!(params.pon_tbi) && params.pon && ('tnscope' in tools || 'mutect2' in tools)) {
        pon_id = pon.map{ it -> [[id:it[0].baseName], it] }
        TABIX_PON(pon_id)
        result_pon_tbi = TABIX_PON.out.tbi.map{ meta, tbi -> [pon, tbi] }
        ch_tabix_versions = ch_tabix_versions.mix(TABIX_PON.out.versions).first()
    }

    ch_versions = ch_versions.mix(ch_tabix_versions)

    result_msisensorpro_scan = Channel.empty()
    if ('msisensorpro' in tools) {
        MSISENSORPRO_SCAN(fasta)
        result_msisensorpro_scan = MSISENSORPRO_SCAN.out.list
        ch_versions = ch_versions.mix(MSISENSORPRO_SCAN.out.versions)
    }

    result_intervals = Channel.empty()
    if (params.no_intervals) {
        file("${params.outdir}/no_intervals.bed").text = "no_intervals\n"
        result_intervals = Channel.fromPath(file("${params.outdir}/no_intervals.bed"))
    } else if (!('annotate' in step) && !('controlfreec' in step)) {
        if (!params.intervals) result_intervals = CREATE_INTERVALS_BED(BUILD_INTERVALS(result_fai))
        else                   result_intervals = CREATE_INTERVALS_BED(file(params.intervals))
    }

    if (!params.no_intervals) {
        result_intervals = result_intervals.flatten()
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

    emit:
        bwa                   = result_bwa
        dbsnp_tbi             = result_dbsnp_tbi
        dict                  = result_dict
        fai                   = result_fai
        germline_resource_tbi = result_germline_resource_tbi
        intervals             = result_intervals
        known_indels_tbi      = result_known_indels_tbi.collect()
        msisensorpro_scan     = result_msisensorpro_scan
        pon_tbi               = result_pon_tbi
        target_bed_gz_tbi     = result_target_bed
        versions              = ch_versions.ifEmpty(null) // channel: [ versions.yml ]
}
