//
// PREPARE VARIANTANNOTATION
//

include { ENSEMBLVEP_DOWNLOAD                       } from '../../../modules/nf-core/ensemblvep/download'
include { SNPEFF_DOWNLOAD                           } from '../../../modules/nf-core/snpeff/download'
include { TABIX_TABIX as TABIX_BCFTOOLS_ANNOTATIONS } from '../../../modules/nf-core/tabix/tabix'

workflow PREPARE_VARIANTANNOTATION {
    take:
    fasta
    bcf_ann_enabled
    bcftools_annotations_params
    bcftools_annotations_tbi_params
    bcftools_header_lines_params
    download_cache
    help_message
    snpeff_enabled
    snpeff_cache_params
    snpeff_genome
    snpeff_db
    vep_enabled
    vep_cache_params
    vep_species
    vep_cache_version
    vep_include_fasta
    vep_genome
    vep_custom_args
    dbnsfp
    dbnsfp_tbi
    spliceai_snv
    spliceai_snv_tbi
    spliceai_indel
    spliceai_indel_tbi

    main:
    bcftools_annotations     = Channel.empty()
    bcftools_annotations_tbi = Channel.empty()
    bcftools_header_lines    = Channel.empty()
    snpeff_cache             = Channel.empty()
    vep_cache                = Channel.empty()
    vep_extra_files          = []
    vep_fasta                = Channel.empty()
    versions                 = Channel.empty()

    if (bcf_ann_enabled) {
        bcftools_annotations  = bcftools_annotations_params  ? Channel.fromPath(bcftools_annotations_params).collect()  : Channel.empty()
        bcftools_header_lines = bcftools_header_lines_params ? Channel.fromPath(bcftools_header_lines_params).collect() : Channel.empty()

        if (bcftools_annotations_tbi_params) {
            bcftools_annotations_tbi = Channel.fromPath(bcftools_annotations_tbi_params).collect()
        } else {
            TABIX_BCFTOOLS_ANNOTATIONS(bcftools_annotations.flatten().map{ it -> [ [ id:it.baseName ], it ] })

            bcftools_annotations_tbi = TABIX_BCFTOOLS_ANNOTATIONS.out.tbi.map{ meta, tbi -> [tbi] }.collect()
            versions = versions.mix(TABIX_BCFTOOLS_ANNOTATIONS.out.versions)
        }
    }

    if (snpeff_enabled) {
        if (download_cache) {
            snpeff_info = Channel.of([ [ id:"${snpeff_genome}.${snpeff_db}" ], snpeff_genome, snpeff_db ])
            SNPEFF_DOWNLOAD(snpeff_info)
            snpeff_cache = SNPEFF_DOWNLOAD.out.cache.collect()
            versions     = versions.mix(SNPEFF_DOWNLOAD.out.versions)
        } else {
            snpeff_annotation_cache_key = (snpeff_cache_params == "s3://annotation-cache/snpeff_cache/") ? "${snpeff_genome}.${snpeff_db}/" : ""
            snpeff_cache_dir =  "${snpeff_annotation_cache_key}${snpeff_genome}.${snpeff_db}"
            snpeff_cache_path_full = file("$snpeff_cache_params/$snpeff_cache_dir", type: 'dir')
            if ( !snpeff_cache_path_full.exists() || !snpeff_cache_path_full.isDirectory() ) {
                if (snpeff_cache_params == "s3://annotation-cache/snpeff_cache/") {
                    error("This path is not available within annotation-cache.\nPlease check https://annotation-cache.github.io/ to create a request for it.")
                } else {
                    error("Path provided with SnpEff cache is invalid.\nMake sure there is a directory named ${snpeff_cache_dir} in ${snpeff_cache_params}./n${help_message}")
                }
            }
            snpeff_cache = Channel.fromPath(file("${snpeff_cache}/${snpeff_annotation_cache_key}"), checkIfExists: true).collect()
                .map{ cache -> [ [ id:"${snpeff_genome}.${snpeff_db}" ], cache ] }
        }
    } else snpeff_cache = []

    if (vep_enabled) {
        vep_fasta = vep_include_fasta && fasta ? Channel.fromPath(fasta).map{ fasta -> [ [id:fasta.baseName], fasta ] }.collect() : [ [id: 'null'], [] ]

        if (dbnsfp && dbnsfp_tbi) {
            vep_extra_files = vep_extra_files.mix(Channel.fromPath(dbnsfp, checkIfExists: true))
            vep_extra_files = vep_extra_files.mix(Channel.fromPath(dbnsfp_tbi, checkIfExists: true))
        }

        if (spliceai_snv && spliceai_snv_tbi && spliceai_indel && spliceai_indel_tbi) {
            vep_extra_files = vep_extra_files.mix(Channel.fromPath(spliceai_indel, checkIfExists: true))
            vep_extra_files = vep_extra_files.mix(Channel.fromPath(spliceai_indel_tbi, checkIfExists: true))
            vep_extra_files = vep_extra_files.mix(Channel.fromPath(spliceai_snv, checkIfExists: true))
            vep_extra_files = vep_extra_files.mix(Channel.fromPath(spliceai_snv_tbi, checkIfExists: true))
        }

        if (download_cache) {
            ensemblvep_info = Channel.of([ [ id:"${vep_cache_version}_${vep_genome}" ], vep_genome, vep_species, vep_cache_version ])
            ENSEMBLVEP_DOWNLOAD(ensemblvep_info)
            vep_cache = ENSEMBLVEP_DOWNLOAD.out.cache.map{ meta, cache -> [ cache ] }.collect()
            versions  = versions.mix(ENSEMBLVEP_DOWNLOAD.out.versions)
        } else {
            vep_annotation_cache_key = (vep_cache_params == "s3://annotation-cache/vep_cache/") ? "${vep_cache_version}_${vep_genome}/" : ""
            vep_species_suffix = vep_custom_args.contains("--merged") ? '_merged' : (vep_custom_args.contains("--refseq") ? '_refseq' : '')
            vep_cache_dir = "${vep_annotation_cache_key}${vep_species}${vep_species_suffix}/${vep_cache_version}_${vep_genome}"
            vep_cache_path_full = file("$vep_cache_params/$vep_cache_dir", type: 'dir')
            if ( !vep_cache_path_full.exists() || !vep_cache_path_full.isDirectory() ) {
                if (vep_cache_params == "s3://annotation-cache/vep_cache/") {
                    error("This path is not available within annotation-cache.\nPlease check https://annotation-cache.github.io/ to create a request for it.")
                } else {
                    error("Path provided with VEP cache is invalid.\nMake sure there is a directory named ${vep_cache_dir} in ${vep_cache_params}./n${help_message}")
                }
            }
            vep_cache = Channel.fromPath(file("${vep_cache}/${vep_annotation_cache_key}"), checkIfExists: true).collect()
        }
    } else vep_cache = []

    emit:
    bcftools_annotations
    bcftools_annotations_tbi // path: bcftools_annotations.vcf.gz.tbi
    bcftools_header_lines
    vep_cache                // channel: [ meta, cache ]
    vep_extra_files
    vep_fasta                // channel: [ meta, fasta ]
    snpeff_cache             // channel: [ meta, cache ]
    versions                 // channel: [ versions.yml ]
}
