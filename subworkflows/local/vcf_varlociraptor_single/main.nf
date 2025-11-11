include { BCFTOOLS_CONCAT as CONCAT_CALLED_CHUNKS   } from '../../../modules/nf-core/bcftools/concat'
include { BCFTOOLS_SORT as SORT_CALLED_CHUNKS       } from '../../../modules/nf-core/bcftools/sort'
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
    ch_versions = Channel.empty()

    meta_map = ch_cram.map { meta, _cram, _crai -> meta +  [sex_string: (meta.sex == "XX" ? "female" : "male")] }

    //TODO this seems suspicious but not the cause for the current resume issues as I am only testing with one sample
    FILL_SCENARIO_FILE(
        meta_map.combine(ch_scenario).map{ meta, scenario_file -> [ meta, scenario_file, [], meta ] }
    )
    ch_scenario_file = FILL_SCENARIO_FILE.out.rendered
    ch_versions = ch_versions.mix(FILL_SCENARIO_FILE.out.versions)


    // Estimate alignment properties
    VARLOCIRAPTOR_ESTIMATEALIGNMENTPROPERTIES(
        ch_cram,
        ch_fasta,
        ch_fasta_fai,
    )
    ch_versions = ch_versions.mix(VARLOCIRAPTOR_ESTIMATEALIGNMENTPROPERTIES.out.versions)

    //
    // CHUNK AND PREPROCESS GERMLINE VCF
    //
    RBT_VCFSPLIT(
        ch_vcf,
        val_num_chunks,
    )
    ch_versions = ch_versions.mix(RBT_VCFSPLIT.out.versions)

    ch_chunked_vcfs = RBT_VCFSPLIT.out.bcfchunks
        .transpose()
        .map { meta, vcf_chunked ->
            [   meta + [chunk: vcf_chunked.name.split(/\./)[-2]],
                vcf_chunked
            ]
        }

    // Join each alignment file with its properties
    ch_cram_alignment = ch_cram.join(VARLOCIRAPTOR_ESTIMATEALIGNMENTPROPERTIES.out.alignment_properties_json, failOnMismatch:true, failOnDuplicate:true)
                                .map{meta, cram, crai, json -> [meta.id, meta, cram, crai, json]}

    // Now combine the each chunked VCFs with the alignment data
    ch_input_preprocess_chunked = ch_chunked_vcfs
        .map { meta, vcf -> [meta.id, meta, vcf] }
        .combine(ch_cram_alignment, by: 0)
        .map { _id, meta_vcf, vcf, meta_cram, cram, crai, alignment_json ->
            [ meta_cram + [
                variantcaller: meta_vcf.variantcaller,
                postprocess: 'varlociraptor',
                chunk: meta_vcf.chunk],
                cram, crai, vcf, alignment_json
            ]
        }

    ch_input_preprocess_chunked.dump(pretty:true, tag: "cram_alignment")

    VARLOCIRAPTOR_PREPROCESS(
        ch_input_preprocess_chunked,
        ch_fasta,
        ch_fasta_fai,
    )
    ch_versions = ch_versions.mix(VARLOCIRAPTOR_PREPROCESS.out.versions)

    //
    // CALL VARIANTS WITH VARLOCIRAPTOR
    //
    VARLOCIRAPTOR_CALLVARIANTS(
        VARLOCIRAPTOR_PREPROCESS.out.bcf,
        ch_scenario_file.map { it -> it[1] }.collect(),
        val_sampletype,
    )
    ch_versions = ch_versions.mix(VARLOCIRAPTOR_CALLVARIANTS.out.versions)

    //
    // SORT AND MERGE CALLED VARIANTS
    //
    SORT_CALLED_CHUNKS(
        VARLOCIRAPTOR_CALLVARIANTS.out.bcf
    )
    ch_versions = ch_versions.mix(SORT_CALLED_CHUNKS.out.versions)

    ch_sort_called_chunks_vcf = SORT_CALLED_CHUNKS.out.vcf
            .branch {
                single: val_num_chunks <= 1
                multiple: val_num_chunks > 1
            }

    ch_sort_called_chunks_tbi = SORT_CALLED_CHUNKS.out.tbi
            .branch {
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

    emit:
    vcf      = ch_sort_called_chunks_vcf.single.mix(CONCAT_CALLED_CHUNKS.out.vcf)
    tbi      = ch_sort_called_chunks_tbi.single.mix(CONCAT_CALLED_CHUNKS.out.tbi)
    versions = ch_versions
}
