include { VARLOCIRAPTOR_CALLVARIANTS                                              } from '../../../modules/nf-core/varlociraptor/callvariants/main'
include { VARLOCIRAPTOR_PREPROCESS as PREPROCESS_TUMOR                            } from '../../../modules/nf-core/varlociraptor/preprocess/main' 
include { VARLOCIRAPTOR_PREPROCESS as PREPROCESS_NORMAL                           } from '../../../modules/nf-core/varlociraptor/preprocess/main'  
include { VARLOCIRAPTOR_ESTIMATEALIGNMENTPROPERTIES as ALIGNMENTPROPERTIES_TUMOR  } from '../../../modules/nf-core/varlociraptor/estimatealignmentproperties/main' 
include { VARLOCIRAPTOR_ESTIMATEALIGNMENTPROPERTIES as ALIGNMENTPROPERTIES_NORMAL } from '../../../modules/nf-core/varlociraptor/estimatealignmentproperties/main' 
include { YTE as FILL_SCENARIO_FILE                                               } from '../../../modules/nf-core/yte/main'
include { RBT_VCFSPLIT as VCFSPLIT_TUMOR                                          } from '../../../modules/nf-core/rbt/vcfsplit/main'
include { RBT_VCFSPLIT as VCFSPLIT_NORMAL                                         } from '../../../modules/nf-core/rbt/vcfsplit/main'
include { BCFTOOLS_SORT as SORT_CALLED_CHUNKS                                     } from '../../../modules/nf-core/bcftools/sort/main'
include { BCFTOOLS_CONCAT as MERGE_CALLED_CHUNKS                                  } from '../../../modules/nf-core/bcftools/concat/main'

workflow VCF_VARLOCIRAPTOR_SOMATIC {

    take:
    ch_cram                          // channel: [mandatory] cram
    ch_fasta                         // channel: [mandatory] fasta
    ch_fasta_fai                     // channel: [mandatory] fasta_fai
    ch_vcf

    main:
    ch_versions = Channel.empty()

    ch_scenario = Channel.fromPath("$projectDir/assets/varlociraptor_somatic_with_priors.yte.yaml", checkIfExists: true)
                    .map{ it -> [ [id:"tumor_normal_scenario"], it ] }

    meta_map = ch_cram.map{ meta, _normal_cram, _normal_crai, _tumor_cram, _tumor_crai -> meta + [sex_string: (meta.sex == "XX" ? "female" : "male") ] }

    // TODO: check if this file is created for each sample, not just once
    FILL_SCENARIO_FILE(
        ch_scenario,
        [],
        meta_map
    )

    ch_scenario_file = FILL_SCENARIO_FILE.out.rendered

    ch_versions = ch_versions.mix(FILL_SCENARIO_FILE.out.versions)

    // TODO: the IDs need to be specific for tumor and normal
    cram_normal = ch_cram.map{ meta, normal_cram, normal_crai, _tumor_cram, _tumor_crai -> [ meta, normal_cram, normal_crai ] }
    cram_tumor  = ch_cram.map{ meta, _normal_cram, _normal_crai, tumor_cram, tumor_crai -> [ meta, tumor_cram, tumor_crai ] }

    // Estimate alignment properties
    ALIGNMENTPROPERTIES_TUMOR(
        cram_normal,
        ch_fasta,
        ch_fasta_fai
    )

    ALIGNMENTPROPERTIES_NORMAL(
        cram_tumor,
        ch_fasta,
        ch_fasta_fai
    )

    ch_versions = ch_versions.mix(ALIGNMENTPROPERTIES_TUMOR.out.versions)
    ch_versions = ch_versions.mix(ALIGNMENTPROPERTIES_NORMAL.out.versions)

    // CHUNK AND PREPROCESS TUMOR VCF
    VCFSPLIT_TUMOR(
        ch_vcf,
        15
    )

    ch_versions = ch_versions.mix(VCFSPLIT_TUMOR.out.versions)

    ch_chunked_tumor_vcfs = VCFSPLIT_TUMOR.out.bcfchunks
        .transpose(by:1)
        .map { meta, vcf_chunked ->
            def new_meta = meta + [chunk:vcf_chunked.name.split(/\./)[-2]]
            [ new_meta, vcf_chunked ]
        }


    // Create a base channel for CRAM data that will be replicated for each chunk
    ch_cram_tumor_base = cram_tumor
        .map{ meta, tumor_cram, tumor_crai -> [ meta.id, meta, tumor_cram, tumor_crai ] }

    // Create a base channel for alignment properties
    ch_alignment_tumor_base = ALIGNMENTPROPERTIES_TUMOR.out.alignment_properties_json
        .map{ meta, alignment_json -> [ meta.id, meta, alignment_json ] }

    // Now combine the chunked VCFs with the base data
    ch_input_tumor_preprocess_chunked = ch_chunked_tumor_vcfs
        .map{ meta, vcf -> [ meta.id, meta, vcf ] }
        .combine(ch_cram_tumor_base, by: 0)  // Combine by sample ID
        .combine(ch_alignment_tumor_base, by: 0)  // Combine by sample ID
        .map{ sample_id, meta_vcf, vcf, meta_cram, tumor_cram, tumor_crai, meta_alignment, alignment_json ->
            def new_meta = meta_cram + [
                variantcaller: meta_vcf.variantcaller, 
                postprocess: 'varlociraptor',
                chunk: meta_vcf.chunk
            ]
            [ new_meta, tumor_cram, tumor_crai, vcf, alignment_json ]
        }

    PREPROCESS_TUMOR(
        ch_input_tumor_preprocess_chunked,
        ch_fasta,
        ch_fasta_fai
    )

    ch_versions = ch_versions.mix(PREPROCESS_TUMOR.out.versions)

    // CHUNK AND PREPROCESS NORMAL VCF
    VCFSPLIT_NORMAL(
        ch_vcf,
        15
    )

    ch_versions = ch_versions.mix(VCFSPLIT_NORMAL.out.versions)

    ch_chunked_normal_vcfs = VCFSPLIT_NORMAL.out.bcfchunks
        .transpose(by:1)
        .map { meta, vcf_chunked ->
            def new_meta = meta + [chunk:vcf_chunked.name.split(/\./)[-2]]
            [ new_meta, vcf_chunked ]
        }



    // Create a base channel for CRAM data that will be replicated for each chunk
    ch_cram_base = cram_normal
        .map{ meta, normal_cram, normal_crai -> [ meta.id, meta, normal_cram, normal_crai ] }

    // Create a base channel for alignment properties
    ch_alignment_base = ALIGNMENTPROPERTIES_NORMAL.out.alignment_properties_json
        .map{ meta, alignment_json -> [ meta.id, meta, alignment_json ] }

    // Now combine the chunked VCFs with the base data
    ch_input_normal_preprocess_chunked = ch_chunked_normal_vcfs
        .map{ meta, vcf -> [ meta.id, meta, vcf ] }
        .combine(ch_cram_base, by: 0)  // Combine by sample ID
        .combine(ch_alignment_base, by: 0)  // Combine by sample ID
        .map{ sample_id, meta_vcf, vcf, meta_cram, normal_cram, normal_crai, meta_alignment, alignment_json ->
            def new_meta = meta_cram + [
                variantcaller: meta_vcf.variantcaller, 
                postprocess: 'varlociraptor',
                chunk: meta_vcf.chunk
            ]
            [ new_meta, normal_cram, normal_crai, vcf, alignment_json ]
        }

    PREPROCESS_NORMAL(
        ch_input_normal_preprocess_chunked,
        ch_fasta,
        ch_fasta_fai
    )

    ch_versions = ch_versions.mix(PREPROCESS_NORMAL.out.versions)

    ch_vcf_for_callvariants = PREPROCESS_NORMAL.out.bcf
            .map{ meta, normal_bcf -> [ meta.id + meta.chunk + meta.variantcaller, meta, normal_bcf]}
        .join(PREPROCESS_TUMOR.out.bcf
            .map{ meta, tumor_bcf -> [ meta.id + meta.chunk + meta.variantcaller, meta, tumor_bcf]}
        ).map{ _id, meta_normal, normal_bcf, _meta_tumor, tumor_bcf ->
            [ meta_normal, [ normal_bcf, tumor_bcf ] ]
        }
        
    VARLOCIRAPTOR_CALLVARIANTS(
        ch_vcf_for_callvariants,
        ch_scenario_file.map{ it -> it[1] }.first(), // scenario file
        Channel.value(["normal", "tumor"]),
    )

    ch_versions = ch_versions.mix(VARLOCIRAPTOR_CALLVARIANTS.out.versions)

    SORT_CALLED_CHUNKS(
        VARLOCIRAPTOR_CALLVARIANTS.out.bcf
    )

    ch_versions = ch_versions.mix(SORT_CALLED_CHUNKS.out.versions)

    ch_vcf_tbi_chunks = SORT_CALLED_CHUNKS.out.vcf
        .map{ meta, vcf -> [[meta.id, meta.variantcaller], meta, vcf] }
        .collate(15)
        .filter{ chunk_list -> 
            //println "VCF Collate: Got ${chunk_list.size()} chunks for ${chunk_list[0][0]}"
            chunk_list.size() == 15 
        }
        .map{ chunk_list ->
            def key = chunk_list[0][0]
            def meta_list = chunk_list.collect{ it[1] }
            def vcf_list = chunk_list.collect{ it[2] }
            [key, meta_list, vcf_list]
        }
        .join(
            SORT_CALLED_CHUNKS.out.tbi
                .map{ meta, tbi -> [[meta.id, meta.variantcaller], meta, tbi] }
                .collate(15)
                .filter{ chunk_list -> 
                    //println "TBI Collate: Got ${chunk_list.size()} chunks for ${chunk_list[0][0]}"
                    chunk_list.size() == 15 
                }
                .map{ chunk_list ->
                    def key = chunk_list[0][0]
                    def meta_list = chunk_list.collect{ it[1] }
                    def tbi_list = chunk_list.collect{ it[2] }
                    [key, meta_list, tbi_list]
                }
        )
        .map{ key, meta_list_vcf, vcf_list, meta_list_tbi, tbi_list ->
            def sample_meta = meta_list_vcf[0].clone()
            sample_meta.remove('chunk')
            [sample_meta, vcf_list, tbi_list]
        }

    MERGE_CALLED_CHUNKS( ch_vcf_tbi_chunks )
    
    ch_final_vcf = MERGE_CALLED_CHUNKS.out.vcf.map{ meta, vcf -> [meta.id + meta.variantcaller, meta, vcf] }
        .join(
            MERGE_CALLED_CHUNKS.out.tbi.map{ meta, tbi -> [meta.id + meta.variantcaller, meta, tbi] }
        )
        .map{ _id, meta_vcf, vcf, _meta_tbi, tbi -> [meta_vcf, vcf, tbi] }
    
    ch_versions = ch_versions.mix(MERGE_CALLED_CHUNKS.out.versions)

    emit:
    vcf      = ch_final_vcf
    versions = ch_versions                     // channel: [ versions.yml ]
}