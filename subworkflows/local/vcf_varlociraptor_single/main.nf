include { BCFTOOLS_CONCAT as CONCAT_CALLED_CHUNKS   } from '../../../modules/nf-core/bcftools/concat'
include { BCFTOOLS_SORT as SORT_CALLED_CHUNKS       } from '../../../modules/nf-core/bcftools/sort'
include { BCFTOOLS_SORT as SORT_FINAL_VCF           } from '../../../modules/nf-core/bcftools/sort'
include { RBT_VCFSPLIT                              } from '../../../modules/nf-core/rbt/vcfsplit'
include { VARLOCIRAPTOR_CALLVARIANTS                } from '../../../modules/nf-core/varlociraptor/callvariants'
include { VARLOCIRAPTOR_ESTIMATEALIGNMENTPROPERTIES } from '../../../modules/nf-core/varlociraptor/estimatealignmentproperties'
include { VARLOCIRAPTOR_PREPROCESS                  } from '../../../modules/nf-core/varlociraptor/preprocess'
include { YTE as FILL_SCENARIO_FILE                 } from '../../../modules/nf-core/yte'

workflow VCF_VARLOCIRAPTOR_SINGLE {
    take:
    ch_cram
    ch_fasta
    ch_fasta_fai
    ch_scenario
    ch_vcf
    val_num_chunks
    val_sampletype

    main:
    ch_versions = channel.empty()

    meta_map = ch_cram.map { meta, _cram, _crai -> meta + [sex_string: (meta.sex == "XX" ? "female" : "male")] }

    FILL_SCENARIO_FILE(
        meta_map.combine(ch_scenario).map { meta, scenario_file -> [meta, scenario_file, [], meta] }
    )
    ch_scenario_file = FILL_SCENARIO_FILE.out.rendered

    // Estimate alignment properties
    VARLOCIRAPTOR_ESTIMATEALIGNMENTPROPERTIES(
        ch_cram.combine(ch_fasta).combine(ch_fasta_fai).map { meta_cram, cram, crai, _meta_fasta, fasta, _meta_fai, fai ->
            [meta_cram, cram, crai, fasta, fai]
        }
    )

    //
    // CHUNK AND PREPROCESS GERMLINE VCF
    //
    RBT_VCFSPLIT(
        ch_vcf,
        val_num_chunks,
    )

    ch_chunked_vcfs = RBT_VCFSPLIT.out.bcfchunks
        .transpose()
        .map { meta, vcf_chunked ->
            [
                meta + [chunk: vcf_chunked.name.split(/\./)[-2]],
                vcf_chunked,
            ]
        }

    // Join each alignment file with its properties
    ch_cram_alignment = ch_cram
        .join(VARLOCIRAPTOR_ESTIMATEALIGNMENTPROPERTIES.out.alignment_properties_json, failOnMismatch: true, failOnDuplicate: true)
        .map { meta, cram, crai, json -> [meta.id, meta, cram, crai, json] }

    // Now combine the each chunked VCFs with the alignment data
    ch_input_preprocess_chunked = ch_chunked_vcfs
        .map { meta, vcf -> [meta.id, meta, vcf] }
        .combine(ch_cram_alignment, by: 0)
        .combine(ch_fasta)
        .combine(ch_fasta_fai)
        .map { _id, meta_vcf, vcf, meta_cram, cram, crai, alignment_json, _meta_fasta, fasta, _meta_fai, fasta_fai ->
            [
                meta_cram + [
                    variantcaller: meta_vcf.variantcaller,
                    postprocess: 'varlociraptor',
                    chunk: meta_vcf.chunk,
                ],
                cram,
                crai,
                vcf,
                alignment_json,
                fasta,
                fasta_fai,
            ]
        }

    VARLOCIRAPTOR_PREPROCESS(
        ch_input_preprocess_chunked
    )

    //
    // CALL VARIANTS WITH VARLOCIRAPTOR
    //
    ch_vcfs_for_callvariants = VARLOCIRAPTOR_PREPROCESS.out.bcf
        .map { meta, bcf ->
            [meta.id, meta, bcf]
        }
        .combine(
            ch_scenario_file.map { meta, scenario_file -> [meta.id, scenario_file] },
            by: 0
        )
        .map { _id, meta_normal, normal_bcf, scenario_file ->
            [meta_normal, [normal_bcf], scenario_file, val_sampletype]
        }

    VARLOCIRAPTOR_CALLVARIANTS(
        ch_vcfs_for_callvariants
    )

    //
    // SORT AND MERGE CALLED VARIANTS
    //
    SORT_CALLED_CHUNKS(
        VARLOCIRAPTOR_CALLVARIANTS.out.bcf
    )
    ch_versions = ch_versions.mix(SORT_CALLED_CHUNKS.out.versions)

    ch_sort_called_chunks_vcf = SORT_CALLED_CHUNKS.out.vcf.branch {
        single: val_num_chunks <= 1
        multiple: val_num_chunks > 1
    }

    ch_sort_called_chunks_tbi = SORT_CALLED_CHUNKS.out.tbi.branch {
        single: val_num_chunks <= 1
        multiple: val_num_chunks > 1
    }

    ch_vcf_tbi_chunks = ch_sort_called_chunks_vcf.multiple
        .join(ch_sort_called_chunks_tbi.multiple, failOnMismatch: true, failOnDuplicate: true)
        .map { meta, vcf, tbi ->
            [meta - meta.subMap("chunk"), vcf, tbi]
        }
        .groupTuple(size: val_num_chunks)

    CONCAT_CALLED_CHUNKS(ch_vcf_tbi_chunks)

    ch_versions = ch_versions.mix(CONCAT_CALLED_CHUNKS.out.versions)

    ch_final_vcf = ch_sort_called_chunks_vcf.single.mix(CONCAT_CALLED_CHUNKS.out.vcf)

    SORT_FINAL_VCF(ch_final_vcf)

    ch_versions = ch_versions.mix(SORT_FINAL_VCF.out.versions)

    emit:
    vcf      = SORT_FINAL_VCF.out.vcf
    tbi      = SORT_FINAL_VCF.out.tbi
    versions = ch_versions
}
