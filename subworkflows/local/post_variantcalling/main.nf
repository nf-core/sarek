//
// POST VARIANT CALLING: processes run on variantcalled but not annotated VCFs
//

include { BCFTOOLS_CONCAT as CONCAT_GERMLINE_STRELKA               } from '../../../modules/nf-core/bcftools/concat'
include { BCFTOOLS_CONCAT as CONCAT_SOMATIC_STRELKA                } from '../../../modules/nf-core/bcftools/concat'
include { CONCATENATE_GERMLINE_VCFS                                } from '../vcf_concatenate_germline'
include { NORMALIZE_VCFS                                           } from '../vcf_normalization'
include { TABIX_TABIX as TABIX_GERMLINE                            } from '../../../modules/nf-core/tabix/tabix'
include { TABIX_TABIX as TABIX_SOMATIC                             } from '../../../modules/nf-core/tabix/tabix'
include { VCF_VARLOCIRAPTOR_SINGLE as VCF_VARLOCIRAPTOR_GERMLINE   } from '../vcf_varlociraptor_single'
include { VCF_VARLOCIRAPTOR_SOMATIC                                } from '../vcf_varlociraptor_somatic'
include { VCF_VARLOCIRAPTOR_SINGLE as VCF_VARLOCIRAPTOR_TUMOR_ONLY } from '../vcf_varlociraptor_single'

workflow POST_VARIANTCALLING {
    take:
    tools
    cram_germline
    germline_vcfs
    cram_tumor_only
    tumor_only_vcfs
    cram_somatic
    somatic_vcfs
    fasta
    fai
    concatenate_vcfs
    normalize_vcfs
    varlociraptor_chunk_size // integer: [mandatory] [default: 15] number of chunks to split BCF files when preprocessing and calling variants

    main:
    versions = Channel.empty()
    vcfs = Channel.empty()

    if (concatenate_vcfs) {
        CONCATENATE_GERMLINE_VCFS(germline_vcfs)

        vcfs = vcfs.mix(CONCATENATE_GERMLINE_VCFS.out.vcfs)
        versions = versions.mix(CONCATENATE_GERMLINE_VCFS.out.versions)
    }

    if (normalize_vcfs) {
        NORMALIZE_VCFS(germline_vcfs, tumor_only_vcfs, somatic_vcfs, fasta)

        vcfs = vcfs.mix(NORMALIZE_VCFS.out.vcfs)

        versions = versions.mix(NORMALIZE_VCFS.out.versions)
    }

    //
    // VARLOCIRAPTOR
    //
    if (tools.split(',').contains('varlociraptor')) {

        TABIX_GERMLINE(germline_vcfs)
        TABIX_SOMATIC(somatic_vcfs)

        germline_vcf_tbi = germline_vcfs.join(TABIX_GERMLINE.out.tbi, by: [0])
        somatic_vcf_tbi = somatic_vcfs.join(TABIX_SOMATIC.out.tbi, by: [0])

        somatic_vcfs.dump(tag: "somatic_vcfs")
        somatic_vcf_tbi.dump(tag: "somatic_vcf_tbi")
        // CONCAT SNV / INDEL VCFs COMING FROM STRELKA
        somatic_branched = somatic_vcf_tbi.branch {
            strelka: it[0].variantcaller == 'strelka'
            other:   it[0].variantcaller != 'strelka'
        }
        germline_branched = germline_vcf_tbi.branch {
            strelka: it[0].variantcaller == 'strelka'
            other:   it[0].variantcaller != 'strelka'
        }
        CONCAT_SOMATIC_STRELKA(somatic_branched.strelka.map { meta, vcf, tbi -> [meta, [vcf], [tbi]] })
        CONCAT_GERMLINE_STRELKA(germline_branched.strelka.map { meta, vcf, tbi -> [meta, [vcf], [tbi]] })

        // GERMLINE
        varlociraptor_scenario_germline = params.varlociraptor_scenario_germline
            ? Channel.fromPath(params.varlociraptor_scenario_germline).map { it -> [[id: it.baseName - '.yte'], it] }.collect()
            : Channel.fromPath("${projectDir}/assets/varlociraptor_germline.yte.yaml").collect()
        VCF_VARLOCIRAPTOR_GERMLINE(cram_germline, fasta, fai, varlociraptor_scenario_germline, germline_vcfs, varlociraptor_chunk_size, 'normal')
        vcfs = vcfs.mix(VCF_VARLOCIRAPTOR_GERMLINE.out.vcf)
        versions = versions.mix(VCF_VARLOCIRAPTOR_GERMLINE.out.versions)

        // SOMATIC
        varlociraptor_scenario_somatic = params.varlociraptor_scenario_somatic
            ? Channel.fromPath(params.varlociraptor_scenario_somatic).map { it -> [[id: it.baseName - '.yte'], it] }.collect()
            : Channel.fromPath("${projectDir}/assets/varlociraptor_somatic.yte.yaml").collect()
        VCF_VARLOCIRAPTOR_SOMATIC(cram_somatic, fasta, fai, varlociraptor_scenario_somatic, somatic_vcf_tbi, germline_vcf_tbi, varlociraptor_chunk_size)
        vcfs = vcfs.mix(VCF_VARLOCIRAPTOR_SOMATIC.out.vcf)
        versions = versions.mix(VCF_VARLOCIRAPTOR_SOMATIC.out.versions)

        // TUMOR ONLY
        varlociraptor_scenario_tumor_only = params.varlociraptor_scenario_tumor_only
            ? Channel.fromPath(params.varlociraptor_scenario_tumor_only).map { it -> [[id: it.baseName - '.yte'], it] }.collect()
            : Channel.fromPath("${projectDir}/assets/varlociraptor_tumor_only.yte.yaml").collect()
        VCF_VARLOCIRAPTOR_TUMOR_ONLY(cram_tumor_only, fasta, fai, varlociraptor_scenario_tumor_only, tumor_only_vcfs, varlociraptor_chunk_size, 'tumor')
        vcfs = vcfs.mix(VCF_VARLOCIRAPTOR_TUMOR_ONLY.out.vcf)
        versions = versions.mix(VCF_VARLOCIRAPTOR_TUMOR_ONLY.out.versions)
    }

    emit:
    vcfs     // post processed vcfs
    versions // channel: [ versions.yml ]
}
