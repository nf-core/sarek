include { BCFTOOLS_CONCAT as CONCAT_CALLED_CHUNKS                                 } from '../../../modules/nf-core/bcftools/concat'
include { BCFTOOLS_CONCAT as CONCAT_SOMATIC_STRELKA                               } from '../../../modules/nf-core/bcftools/concat'
include { BCFTOOLS_MERGE as MERGE_GERMLINE_SOMATIC_VCFS                           } from '../../../modules/nf-core/bcftools/merge'
include { BCFTOOLS_SORT as SORT_CALLED_CHUNKS                                     } from '../../../modules/nf-core/bcftools/sort'
include { BCFTOOLS_SORT as SORT_FINAL_VCF                                         } from '../../../modules/nf-core/bcftools/sort'
include { TABIX_TABIX as TABIX_GERMLINE                                           } from '../../../modules/nf-core/tabix/tabix'
include { TABIX_TABIX as TABIX_SOMATIC                                            } from '../../../modules/nf-core/tabix/tabix'
include { RBT_VCFSPLIT                                                            } from '../../../modules/nf-core/rbt/vcfsplit'
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
    ch_versions = channel.empty()

    meta_map = ch_cram.map { meta, _normal_cram, _normal_crai, _tumor_cram, _tumor_crai ->
        meta + [sex_string: (meta.sex == "XX" ? "female" : "male")]
    }

    FILL_SCENARIO_FILE(
        meta_map.combine(ch_scenario).map { meta, scenario_file -> [meta, scenario_file, [], meta] }
    )
    ch_scenario_file = FILL_SCENARIO_FILE.out.rendered

    cram_normal = ch_cram.map { meta, normal_cram, normal_crai, _tumor_cram, _tumor_crai -> [meta + [match_id: meta.normal_id], normal_cram, normal_crai] }
    cram_tumor = ch_cram.map { meta, _normal_cram, _normal_crai, tumor_cram, tumor_crai -> [meta + [match_id: meta.normal_id], tumor_cram, tumor_crai] }

    // Estimate alignment properties
    ALIGNMENTPROPERTIES_TUMOR(
        cram_tumor.combine(ch_fasta).combine(ch_fasta_fai).map { meta_cram, cram, crai, _meta_fasta, fasta, _meta_fai, fai ->
            [meta_cram, cram, crai, fasta, fai]
        }
    )

    ALIGNMENTPROPERTIES_NORMAL(
        cram_normal.combine(ch_fasta).combine(ch_fasta_fai).map { meta_cram, cram, crai, _meta_fasta, fasta, _meta_fai, fai ->
            [meta_cram, cram, crai, fasta, fai]
        }
    )

    //
    // CONCAT SNV AND INDEL VCFS FOR STRELKA
    //
    TABIX_SOMATIC(ch_somatic_vcf)
    ch_versions = ch_versions.mix(TABIX_SOMATIC.out.versions)
    ch_somatic_vcf_tbi = ch_somatic_vcf.join(TABIX_SOMATIC.out.tbi, by: [0])

    // CONCAT SNV / INDEL VCFs COMING FROM STRELKA
    ch_somatic_branched = ch_somatic_vcf_tbi.branch { items ->
        strelka: items[0].variantcaller == 'strelka'
        other: items[0].variantcaller != 'strelka'
    }

    // Group somatic strelka SNVs and INDELs by sample for concatenation
    ch_somatic_strelka_grouped = ch_somatic_branched.strelka
        .map { meta, vcf, tbi -> [[meta.normal_id, meta.patient], meta, vcf, tbi] }
        .groupTuple(by: 0)
        .map { _key, meta, vcf_list, tbi_list ->
            [meta[0], vcf_list, tbi_list]
        }

    CONCAT_SOMATIC_STRELKA(ch_somatic_strelka_grouped)

    // Use concatenated Strelka VCFs for somatic and germline calling, mix with other variant callers
    ch_somatic_vcf_conc = CONCAT_SOMATIC_STRELKA.out.vcf
        .join(CONCAT_SOMATIC_STRELKA.out.tbi, by: [0])
        .mix(ch_somatic_branched.other)

    //
    // MERGE GERMLINE AND SOMATIC VCFs
    //
    TABIX_GERMLINE(ch_germline_vcf)
    ch_versions = ch_versions.mix(TABIX_GERMLINE.out.versions)
    ch_germline_vcf_tbi = ch_germline_vcf.join(TABIX_GERMLINE.out.tbi, by: [0])

    def somatic_with_key = ch_somatic_vcf_conc.map { meta, vcf, tbi ->
        [[id: meta.normal_id, variantcaller: meta.variantcaller], meta, vcf, tbi]
    }

    def germline_with_key = ch_germline_vcf_tbi.map { meta, vcf, tbi ->
        [[id: meta.id, variantcaller: meta.variantcaller], meta, vcf, tbi]
    }

    def matching_pairs = somatic_with_key.join(germline_with_key, failOnMismatch: false)

    // Branch based on whether a matching germline VCF was found
    def branched = matching_pairs.branch { items ->
        matched: items.size() == 7
        unmatched: items.size() == 4
    }

    MERGE_GERMLINE_SOMATIC_VCFS(
        branched.matched.map { _key, meta_somatic, somatic_vcf, somatic_tbi, _meta_germline, germline_vcf, germline_tbi ->
            [meta_somatic, [somatic_vcf, germline_vcf], [somatic_tbi, germline_tbi]]
        },
        ch_fasta,
        ch_fasta_fai,
        [[], []],
    )

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

    //
    // SPLIT VCF CHUNKS - create chunked VCFs for both tumor and normal preprocessing
    //
    ch_chunked_tumor_vcfs = RBT_VCFSPLIT.out.bcfchunks
        .transpose()
        .map { meta, vcf_chunked ->
            [
                meta + [chunk: vcf_chunked.name.split(/\./)[-2]],
                vcf_chunked,
            ]
        }

    // Create chunked VCFs with normal metadata for normal preprocessing
    ch_chunked_normal_vcfs = RBT_VCFSPLIT.out.bcfchunks
        .transpose()
        .map { meta, vcf_chunked ->
            [
                meta + [match_id: meta.normal_id, chunk: vcf_chunked.name.split(/\./)[-2]],
                vcf_chunked,
            ]
        }

    //
    // PREPROCESS VCF WITH TUMOR CRAM
    //

    // Create base channels for data that will be replicated for each chunk
    ch_cram_tumor = cram_tumor
        .join(ALIGNMENTPROPERTIES_TUMOR.out.alignment_properties_json, failOnMismatch: true, failOnDuplicate: true)
        .map { meta, cram, crai, json -> [meta.id, meta, cram, crai, json] }

    ch_input_tumor_preprocess_chunked = ch_chunked_tumor_vcfs
        .map { meta, vcf -> [meta.id, meta, vcf] }
        .combine(ch_cram_tumor, by: 0)
        .combine(ch_fasta)
        .combine(ch_fasta_fai)
        .map { _id, meta_vcf, vcf, meta_cram, tumor_cram, tumor_crai, alignment_json, _meta_fasta, fasta, _meta_fai, fai ->
            [
                meta_cram + [
                    variantcaller: meta_vcf.variantcaller,
                    postprocess: 'varlociraptor',
                    chunk: meta_vcf.chunk,
                ],
                tumor_cram,
                tumor_crai,
                vcf,
                alignment_json,
                fasta,
                fai,
            ]
        }

    PREPROCESS_TUMOR(
        ch_input_tumor_preprocess_chunked
    )

    //
    // PREPROCESS VCF WITH NORMAL CRAM
    //

    // Create base channels for data that will be replicated for each chunk
    ch_cram_alignment = cram_normal
        .join(ALIGNMENTPROPERTIES_NORMAL.out.alignment_properties_json, failOnMismatch: true, failOnDuplicate: true)
        .map { meta, cram, crai, json -> [meta.match_id, meta, cram, crai, json] }

    ch_input_normal_preprocess_chunked = ch_chunked_normal_vcfs
        .map { meta, vcf -> [meta.match_id, meta, vcf] }
        .combine(ch_cram_alignment, by: 0)
        .combine(ch_fasta)
        .combine(ch_fasta_fai)
        .map { _match_id, meta_vcf, vcf, meta_cram, normal_cram, normal_crai, alignment_json, _meta_fasta, fasta, _meta_fai, fai ->
            [
                meta_vcf + [
                    id: meta_cram.id,
                    postprocess: 'varlociraptor',
                ],
                normal_cram,
                normal_crai,
                vcf,
                alignment_json,
                fasta,
                fai,
            ]
        }

    PREPROCESS_NORMAL(
        ch_input_normal_preprocess_chunked
    )

    //
    // CALL VARIANTS WITH VARLOCIRAPTOR
    //
    ch_normal_for_join = PREPROCESS_NORMAL.out.bcf
        .map { meta, normal_bcf -> [[meta.patient, meta.match_id, meta.chunk, meta.variantcaller], meta, normal_bcf] }
        .dump(tag: 'NORMAL_KEY', pretty: true)

    ch_tumor_for_join = PREPROCESS_TUMOR.out.bcf
        .map { meta, tumor_bcf -> [[meta.patient, meta.normal_id, meta.chunk, meta.variantcaller], meta, tumor_bcf] }
        .dump(tag: 'TUMOR_KEY', pretty: true)

    ch_vcf_for_callvariants = ch_normal_for_join
        .join(
            ch_tumor_for_join,
            by: [0],
            failOnMismatch: true,
            failOnDuplicate: true,
        )
        .combine(ch_scenario_file)
        .map { _id, meta_normal, normal_bcf, _meta_tumor, tumor_bcf, _meta_scenario, scenario_file ->
            [meta_normal, [normal_bcf, tumor_bcf], scenario_file, ["normal", "tumor"]]
        }

    VARLOCIRAPTOR_CALLVARIANTS(
        ch_vcf_for_callvariants
    )

    //
    // SORT AND MERGE CALLED VARIANTS
    //
    SORT_CALLED_CHUNKS(
        VARLOCIRAPTOR_CALLVARIANTS.out.bcf
    )

    ch_sort_called_chunks_vcf = SORT_CALLED_CHUNKS.out.vcf.branch {
        single: val_num_chunks <= 1
        multiple: val_num_chunks > 1
    }

    ch_sort_called_chunks_tbi = SORT_CALLED_CHUNKS.out.tbi.branch {
        single: val_num_chunks <= 1
        multiple: val_num_chunks > 1
    }

    ch_vcf_tbi_chunks = ch_sort_called_chunks_vcf.multiple
        .join(ch_sort_called_chunks_tbi.multiple, by: 0, failOnMismatch: true, failOnDuplicate: true)
        .map { meta, vcf, tbi ->
            [meta - meta.subMap("chunk"), vcf, tbi]
        }
        .groupTuple(size: val_num_chunks)

    CONCAT_CALLED_CHUNKS(ch_vcf_tbi_chunks)

    ch_final_vcf = ch_sort_called_chunks_vcf.single.mix(CONCAT_CALLED_CHUNKS.out.vcf)

    SORT_FINAL_VCF(ch_final_vcf)

    emit:
    vcf      = SORT_FINAL_VCF.out.vcf
    tbi      = SORT_FINAL_VCF.out.tbi
    versions = ch_versions
}
