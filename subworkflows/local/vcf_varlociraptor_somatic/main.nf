include { BCFTOOLS_CONCAT as MERGE_CALLED_CHUNKS                                  } from '../../../modules/nf-core/bcftools/concat'
include { BCFTOOLS_MERGE as MERGE_GERMLINE_SOMATIC_VCFS                           } from '../../../modules/nf-core/bcftools/merge'
include { BCFTOOLS_SORT as SORT_CALLED_CHUNKS                                     } from '../../../modules/nf-core/bcftools/sort'
include { RBT_VCFSPLIT                                                            } from '../../../modules/nf-core/rbt/vcfsplit'
include { TABIX_TABIX as TABIX_GERMLINE                                           } from '../../../modules/nf-core/tabix/tabix'
include { TABIX_TABIX as TABIX_SOMATIC                                            } from '../../../modules/nf-core/tabix/tabix'
include { VARLOCIRAPTOR_CALLVARIANTS                                              } from '../../../modules/nf-core/varlociraptor/callvariants'
include { VARLOCIRAPTOR_ESTIMATEALIGNMENTPROPERTIES as ALIGNMENTPROPERTIES_NORMAL } from '../../../modules/nf-core/varlociraptor/estimatealignmentproperties'
include { VARLOCIRAPTOR_ESTIMATEALIGNMENTPROPERTIES as ALIGNMENTPROPERTIES_TUMOR  } from '../../../modules/nf-core/varlociraptor/estimatealignmentproperties'
include { VARLOCIRAPTOR_PREPROCESS as PREPROCESS_NORMAL                           } from '../../../modules/nf-core/varlociraptor/preprocess'
include { VARLOCIRAPTOR_PREPROCESS as PREPROCESS_TUMOR                            } from '../../../modules/nf-core/varlociraptor/preprocess'
include { YTE as FILL_SCENARIO_FILE                                               } from '../../../modules/nf-core/yte'

workflow VCF_VARLOCIRAPTOR_SOMATIC {
    take:
    ch_cram
    ch_fasta
    ch_fasta_fai
    ch_scenario
    ch_somatic_vcf
    ch_germline_vcf
    val_num_chunks

    main:
    ch_versions = Channel.empty()

    meta_map = ch_cram.map { meta, _normal_cram, _normal_crai, _tumor_cram, _tumor_crai ->
        meta + [sex_string: (meta.sex == "XX" ? "female" : "male")]
    }

    FILL_SCENARIO_FILE(
        meta_map.combine(ch_scenario),
        [],
        meta_map,
    )

    ch_scenario_file = FILL_SCENARIO_FILE.out.rendered

    ch_versions = ch_versions.mix(FILL_SCENARIO_FILE.out.versions)

    cram_normal = ch_cram.map { meta, normal_cram, normal_crai, _tumor_cram, _tumor_crai -> [meta, normal_cram, normal_crai] }
    cram_tumor = ch_cram.map { meta, _normal_cram, _normal_crai, tumor_cram, tumor_crai -> [meta, tumor_cram, tumor_crai] }

    // Estimate alignment properties
    ALIGNMENTPROPERTIES_TUMOR(
        cram_normal,
        ch_fasta,
        ch_fasta_fai,
    )

    ALIGNMENTPROPERTIES_NORMAL(
        cram_tumor,
        ch_fasta,
        ch_fasta_fai,
    )

    ch_versions = ch_versions.mix(ALIGNMENTPROPERTIES_TUMOR.out.versions)
    ch_versions = ch_versions.mix(ALIGNMENTPROPERTIES_NORMAL.out.versions)

    //
    // CONCAT GERMLINE AND SOMATIC VCFs
    //
    TABIX_GERMLINE(ch_germline_vcf)
    TABIX_SOMATIC(ch_somatic_vcf)

    ch_germline_vcf_tbi = ch_germline_vcf.join(TABIX_GERMLINE.out.tbi, by: [0])
    ch_somatic_vcf_tbi = ch_somatic_vcf.join(TABIX_SOMATIC.out.tbi, by: [0])

    ch_versions = ch_versions.mix(TABIX_GERMLINE.out.versions)
    ch_versions = ch_versions.mix(TABIX_SOMATIC.out.versions)

    def somatic_with_key = ch_somatic_vcf_tbi.map { meta, vcf, tbi ->
        [[id: meta.normal_id, variantcaller: meta.variantcaller], meta, vcf, tbi]
    }

    def germline_with_key = ch_germline_vcf_tbi.map { meta, vcf, tbi ->
        [[id: meta.id, variantcaller: meta.variantcaller], meta, vcf, tbi]
    }

    def matching_pairs = somatic_with_key.join(germline_with_key, failOnMismatch: false)

    // Branch based on whether a matching germline VCF was found
    def branched = matching_pairs.branch {
        matched: it.size() == 7
        unmatched: it.size() == 4
    }

    branched.matched.dump(tag: "matched_pairs")
    branched.unmatched.dump(tag: "unmatched_somatic")

    MERGE_GERMLINE_SOMATIC_VCFS(
        branched.matched.map { _key, meta_somatic, somatic_vcf, somatic_tbi, _meta_germline, germline_vcf, germline_tbi ->
            [meta_somatic, [somatic_vcf, germline_vcf], [somatic_tbi, germline_tbi]]
        },
        ch_fasta,
        ch_fasta_fai,
        [[], []],
    )

    ch_versions = ch_versions.mix(MERGE_GERMLINE_SOMATIC_VCFS.out.versions)

    // Combine merged VCFs with unmatched somatic VCFs
    ch_vcf = MERGE_GERMLINE_SOMATIC_VCFS.out.vcf.mix(
        branched.unmatched.map { _key, meta, vcf, _tbi -> [meta, vcf] }
    )

    //
    // CHUNK VCF FILES
    //
    RBT_VCFSPLIT(
        ch_vcf,
        val_num_chunks,
    )

    ch_versions = ch_versions.mix(RBT_VCFSPLIT.out.versions)

    //
    // PREPROCESS VCF WITH TUMOR CRAM
    //
    ch_chunked_tumor_vcfs = RBT_VCFSPLIT.out.bcfchunks
        .transpose(by: 1)
        .map { meta, vcf_chunked ->
            def new_meta = meta + [chunk: vcf_chunked.name.split(/\./)[-2]]
            [new_meta, vcf_chunked]
        }

    // Create base channels for data that will be replicated for each chunk
    ch_cram_tumor_base = cram_tumor.map { meta, tumor_cram, tumor_crai -> [meta.id, meta, tumor_cram, tumor_crai] }

    ch_alignment_tumor_base = ALIGNMENTPROPERTIES_TUMOR.out.alignment_properties_json.map { meta, alignment_json -> [meta.id, meta, alignment_json] }

    ch_input_tumor_preprocess_chunked = ch_chunked_tumor_vcfs
        .map { meta, vcf -> [meta.id, meta, vcf] }
        .combine(ch_cram_tumor_base, by: 0)
        .combine(ch_alignment_tumor_base, by: 0)
        .map { _id, meta_vcf, vcf, meta_cram, tumor_cram, tumor_crai, _meta_alignment, alignment_json ->
            def new_meta = meta_cram + [
                variantcaller: meta_vcf.variantcaller,
                postprocess: 'varlociraptor',
                chunk: meta_vcf.chunk,
            ]
            [new_meta, tumor_cram, tumor_crai, vcf, alignment_json]
        }

    PREPROCESS_TUMOR(
        ch_input_tumor_preprocess_chunked,
        ch_fasta,
        ch_fasta_fai,
    )

    ch_versions = ch_versions.mix(PREPROCESS_TUMOR.out.versions)

    //
    // PREPROCESS VCF WITH NORMAL CRAM
    //
    ch_chunked_normal_vcfs = RBT_VCFSPLIT.out.bcfchunks
        .transpose(by: 1)
        .map { meta, vcf_chunked ->
            def new_meta = meta + [chunk: vcf_chunked.name.split(/\./)[-2]]
            [new_meta, vcf_chunked]
        }

    // Create base channels for data that will be replicated for each chunk
    ch_cram_base = cram_normal.map { meta, normal_cram, normal_crai -> [meta.id, meta, normal_cram, normal_crai] }

    ch_alignment_base = ALIGNMENTPROPERTIES_NORMAL.out.alignment_properties_json.map { meta, alignment_json -> [meta.id, meta, alignment_json] }

    ch_input_normal_preprocess_chunked = ch_chunked_normal_vcfs
        .map { meta, vcf -> [meta.id, meta, vcf] }
        .combine(ch_cram_base, by: 0)
        .combine(ch_alignment_base, by: 0)
        .map { _id, meta_vcf, vcf, meta_cram, normal_cram, normal_crai, _meta_alignment, alignment_json ->
            def new_meta = meta_cram + [
                variantcaller: meta_vcf.variantcaller,
                postprocess: 'varlociraptor',
                chunk: meta_vcf.chunk,
            ]
            [new_meta, normal_cram, normal_crai, vcf, alignment_json]
        }

    PREPROCESS_NORMAL(
        ch_input_normal_preprocess_chunked,
        ch_fasta,
        ch_fasta_fai,
    )

    ch_versions = ch_versions.mix(PREPROCESS_NORMAL.out.versions)

    //
    // CALL VARIANTS WITH VARLOCIRAPTOR
    //
    ch_vcf_for_callvariants = PREPROCESS_NORMAL.out.bcf
        .map { meta, normal_bcf -> [meta.id + meta.chunk + meta.variantcaller, meta, normal_bcf] }
        .join(
            PREPROCESS_TUMOR.out.bcf.map { meta, tumor_bcf -> [meta.id + meta.chunk + meta.variantcaller, meta, tumor_bcf] }
        )
        .map { _id, meta_normal, normal_bcf, _meta_tumor, tumor_bcf ->
            [meta_normal, [normal_bcf, tumor_bcf]]
        }

    VARLOCIRAPTOR_CALLVARIANTS(
        ch_vcf_for_callvariants,
        ch_scenario_file.map { it -> it[1] }.first(),
        Channel.value(["normal", "tumor"]),
    )

    ch_versions = ch_versions.mix(VARLOCIRAPTOR_CALLVARIANTS.out.versions)

    //
    // SORT AND MERGE CALLED VARIANTS
    //
    SORT_CALLED_CHUNKS(VARLOCIRAPTOR_CALLVARIANTS.out.bcf)

    ch_versions = ch_versions.mix(SORT_CALLED_CHUNKS.out.versions)

    ch_vcf_tbi_chunks = SORT_CALLED_CHUNKS.out.vcf
        .join(SORT_CALLED_CHUNKS.out.tbi, failOnMismatch: true, failOnDuplicate: true)
        .map { meta, vcf, tbi ->
            def new_meta = meta - meta.subMap("chunk")
            [new_meta, vcf, tbi]
        }
        .groupTuple(size: val_num_chunks)

    MERGE_CALLED_CHUNKS(ch_vcf_tbi_chunks)

    ch_final_vcf = MERGE_CALLED_CHUNKS.out.vcf
        .map { meta, vcf -> [meta.id + meta.variantcaller, meta, vcf] }
        .join(
            MERGE_CALLED_CHUNKS.out.tbi.map { meta, tbi -> [meta.id + meta.variantcaller, meta, tbi] }
        )
        .map { _id, meta_vcf, vcf, _meta_tbi, tbi -> [meta_vcf, vcf, tbi] }

    ch_versions = ch_versions.mix(MERGE_CALLED_CHUNKS.out.versions)

    emit:
    vcf      = ch_final_vcf
    versions = ch_versions
}
