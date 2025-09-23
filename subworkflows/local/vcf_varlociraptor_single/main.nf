include { BCFTOOLS_CONCAT as MERGE_CALLED_CHUNKS    } from '../../../modules/nf-core/bcftools/concat'
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
    // 'normal' or 'tumor'
    val_sampletype

    main:
    ch_versions = Channel.empty()

    meta_map = ch_cram.map { meta, _cram, _crai -> meta + [sex_string: (meta.sex == "XX" ? "female" : "male")] }

    FILL_SCENARIO_FILE(
        meta_map.combine(ch_scenario),
        [],
        meta_map,
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

    ch_chunked_tumor_vcfs = RBT_VCFSPLIT.out.bcfchunks
        .transpose(by: 1)
        .map { meta, vcf_chunked ->
            def new_meta = meta + [chunk: vcf_chunked.name.split(/\./)[-2]]
            [new_meta, vcf_chunked]
        }

    // Create a base channel for CRAM data that will be replicated for each chunk
    ch_cram_base = ch_cram.map { meta, cram, crai -> [meta.id, meta, cram, crai] }

    // Create a base channel for alignment properties
    ch_alignment_base = VARLOCIRAPTOR_ESTIMATEALIGNMENTPROPERTIES.out.alignment_properties_json.map { meta, alignment_json -> [meta.id, meta, alignment_json] }

    // Now combine the chunked VCFs with the base data
    ch_input_tumor_preprocess_chunked = ch_chunked_tumor_vcfs
        .map { meta, vcf -> [meta.id, meta, vcf] }
        .combine(ch_cram_base, by: 0)
        .combine(ch_alignment_base, by: 0)
        .map { _id, meta_vcf, vcf, meta_cram, cram, crai, _meta_alignment, alignment_json ->
            def new_meta = meta_cram + [
                variantcaller: meta_vcf.variantcaller,
                postprocess: 'varlociraptor',
                chunk: meta_vcf.chunk,
            ]
            [new_meta, cram, crai, vcf, alignment_json]
        }

    VARLOCIRAPTOR_PREPROCESS(
        ch_input_tumor_preprocess_chunked,
        ch_fasta,
        ch_fasta_fai,
    )

    ch_versions = ch_versions.mix(VARLOCIRAPTOR_PREPROCESS.out.versions)

    //
    // CALL VARIANTS WITH VARLOCIRAPTOR
    //
    VARLOCIRAPTOR_CALLVARIANTS(
        VARLOCIRAPTOR_PREPROCESS.out.bcf,
        ch_scenario_file.map { it -> it[1] }.first(),
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
