include { VARLOCIRAPTOR_CALLVARIANTS                                              } from '../../../modules/nf-core/varlociraptor/callvariants/main'
include { VARLOCIRAPTOR_PREPROCESS as PREPROCESS_TUMOR                            } from '../../../modules/nf-core/varlociraptor/preprocess/main' 
include { VARLOCIRAPTOR_PREPROCESS as PREPROCESS_NORMAL                           } from '../../../modules/nf-core/varlociraptor/preprocess/main'  
include { VARLOCIRAPTOR_ESTIMATEALIGNMENTPROPERTIES as ALIGNMENTPROPERTIES_TUMOR  } from '../../../modules/nf-core/varlociraptor/estimatealignmentproperties/main' 
include { VARLOCIRAPTOR_ESTIMATEALIGNMENTPROPERTIES as ALIGNMENTPROPERTIES_NORMAL } from '../../../modules/nf-core/varlociraptor/estimatealignmentproperties/main' 
include { YTE as FILL_SCENARIO_FILE                                               } from '../../../modules/nf-core/yte/main'

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

    ch_input_normal_preprocess = cram_normal
        .map{ meta, normal_cram, normal_crai -> [ [id: meta.id], meta, normal_cram, normal_crai ] }
        .join(
            ch_vcf
                .map{ meta, vcf -> [ [id: meta.id], meta , vcf ] }
        ).join(
            ALIGNMENTPROPERTIES_NORMAL.out.alignment_properties_json
                .map{ meta, alignment_json -> [ [id: meta.id], meta, alignment_json ] }
        ).map{
            _id, meta_cram, normal_cram, normal_crai, meta_vcf, vcf, _meta_json, alignment_json ->
            def new_meta = meta_cram + [variantcaller: meta_vcf.variantcaller, postprocess: 'varlociraptor']
            [ new_meta, normal_cram, normal_crai, vcf, alignment_json ]
        }

    ch_input_tumor_preprocess = cram_tumor
        .map{ meta, tumor_cram, tumor_crai -> [ [id: meta.id], meta, tumor_cram, tumor_crai ] }
        .join(
            ch_vcf
                .map{ meta, vcf -> [ [id: meta.id], meta , vcf ] }
        ).join(
            ALIGNMENTPROPERTIES_NORMAL.out.alignment_properties_json
                .map{ meta, alignment_json -> [ [id: meta.id], meta, alignment_json ] }
        ).map{
            _id, meta_cram, tumor_cram, tumor_crai, meta_vcf, vcf, _meta_json, alignment_json ->
            def new_meta = meta_cram + [variantcaller: meta_vcf.variantcaller, postprocess: 'varlociraptor']
            [ new_meta, tumor_cram, tumor_crai, vcf, alignment_json ]
        }

    // TODO: add chunking here

    PREPROCESS_NORMAL(
        ch_input_normal_preprocess,
        ch_fasta,
        ch_fasta_fai
    )

    PREPROCESS_TUMOR(
        ch_input_tumor_preprocess,
        ch_fasta,
        ch_fasta_fai
    )

    // TODO: add sort & split here

    ch_versions = ch_versions.mix(PREPROCESS_NORMAL.out.versions)
    ch_versions = ch_versions.mix(PREPROCESS_TUMOR.out.versions)

    // scenario_list ("tumor", "normal")
    ch_vcf_for_callvariants = PREPROCESS_NORMAL.out.bcf
        .join(PREPROCESS_TUMOR.out.bcf).map{ meta_normal, normal_bcf, _meta_tumor, tumor_bcf ->
            [ meta_normal, [ normal_bcf, tumor_bcf ] ]
        }

    VARLOCIRAPTOR_CALLVARIANTS(
        ch_vcf_for_callvariants,
        ch_scenario_file,
        ["normal", "tumor"] 
    )

    ch_versions = ch_versions.mix(VARLOCIRAPTOR_CALLVARIANTS.out.versions)

    // TODO: group (patient id) (welche samples werden zusammengefasst?)
    emit:
    vcf      = VARLOCIRAPTOR_CALLVARIANTS.out.bcf
    versions = ch_versions                     // channel: [ versions.yml ]
}