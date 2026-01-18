//
// POST VARIANT CALLING: processes run on variantcalled but not annotated VCFs
//
include { BCFTOOLS_VIEW as FILTER_VCFS                             } from '../../../modules/nf-core/bcftools/view'
include { CONCATENATE_GERMLINE_VCFS                                } from '../vcf_concatenate_germline'
include { CONSENSUS                                                } from '../vcf_consensus'
include { NORMALIZE_VCFS                                           } from '../vcf_normalization'
include { VCF_VARLOCIRAPTOR_SINGLE as VCF_VARLOCIRAPTOR_GERMLINE   } from '../vcf_varlociraptor_single'
include { VCF_VARLOCIRAPTOR_SOMATIC                                } from '../vcf_varlociraptor_somatic'
include { VCF_VARLOCIRAPTOR_SINGLE as VCF_VARLOCIRAPTOR_TUMOR_ONLY } from '../vcf_varlociraptor_single'

workflow POST_VARIANTCALLING {
    take:
    tools
    cram_germline
    germline_vcfs
    germline_tbis
    cram_tumor_only
    tumor_only_vcfs
    tumor_only_tbis
    cram_somatic
    somatic_vcfs
    somatic_tbis
    fasta
    fai
    concatenate_vcfs
    filter_vcfs
    snv_consensus_calling
    normalize_vcfs
    varlociraptor_chunk_size // integer: [mandatory] [default: 15] number of chunks to split BCF files when preprocessing and calling variants
    varlociraptor_scenario_germline
    varlociraptor_scenario_somatic
    varlociraptor_scenario_tumor_only

    main:
    versions = Channel.empty()
    vcfs = Channel.empty()
    tbis = Channel.empty()

    //
    // VARLOCIRAPTOR
    //
    if (tools && tools.split(',').contains('varlociraptor')) {
        // GERMLINE
        VCF_VARLOCIRAPTOR_GERMLINE(cram_germline, fasta, fai, varlociraptor_scenario_germline, germline_vcfs, varlociraptor_chunk_size, 'normal')

        vcfs = vcfs.mix(VCF_VARLOCIRAPTOR_GERMLINE.out.vcf)
        tbis = tbis.mix(VCF_VARLOCIRAPTOR_GERMLINE.out.tbi)
        versions = versions.mix(VCF_VARLOCIRAPTOR_GERMLINE.out.versions)

        // SOMATIC
        VCF_VARLOCIRAPTOR_SOMATIC(cram_somatic, fasta, fai, varlociraptor_scenario_somatic, somatic_vcfs, germline_vcfs, varlociraptor_chunk_size)

        vcfs = vcfs.mix(VCF_VARLOCIRAPTOR_SOMATIC.out.vcf)
        tbis = tbis.mix(VCF_VARLOCIRAPTOR_SOMATIC.out.tbi)
        versions = versions.mix(VCF_VARLOCIRAPTOR_SOMATIC.out.versions)

        // TUMOR ONLY
        VCF_VARLOCIRAPTOR_TUMOR_ONLY(cram_tumor_only, fasta, fai, varlociraptor_scenario_tumor_only, tumor_only_vcfs, varlociraptor_chunk_size, 'tumor')

        vcfs = vcfs.mix(VCF_VARLOCIRAPTOR_TUMOR_ONLY.out.vcf)
        tbis = tbis.mix(VCF_VARLOCIRAPTOR_TUMOR_ONLY.out.tbi)
        versions = versions.mix(VCF_VARLOCIRAPTOR_TUMOR_ONLY.out.versions)

    } else if (filter_vcfs || normalize_vcfs || concatenate_vcfs ) {

        // IMPORTANT: When adding new SNV variant callers to Sarek, add them to this list!
        // This list determines which variant callers are eligible for:
        // - VCF normalization (--normalize_vcfs)
        // - VCF filtering (--filter_vcfs)
        // - Consensus calling (--snv_consensus_calling)
        //
        // To find all variant callers: grep "variantcaller:" subworkflows/local/bam_variant_calling*/main.nf
        //
        // Excluded callers (not eligible for normalization/consensus):
        // - manta, tiddit: structural variant callers (separate workflow)
        // - samtools mpileup produces pileup format for ControlFREEC, not consensus-ready VCFs
        def small_variantcallers = ['bcftools', 'deepvariant', 'freebayes', 'haplotypecaller',
                                    'lofreq', 'muse', 'mutect2', 'sentieon_dnascope',
                                    'sentieon_haplotyper', 'sentieon_tnscope', 'strelka' ]

        def excluded_variantcallers = ['manta', 'tiddit', 'samtools']

        all_vcfs = Channel.empty().mix(germline_vcfs, tumor_only_vcfs, somatic_vcfs)
                                .branch{ meta, vcf ->
                                    small: small_variantcallers.contains(meta.variantcaller)
                                    other: true
                                }

        all_tbis = Channel.empty().mix(germline_tbis, tumor_only_tbis, somatic_tbis)
                                .branch{ meta, tbi ->
                                    small: small_variantcallers.contains(meta.variantcaller)
                                    other: true
                                }

        // Validate that we're not silently excluding unknown variant callers
        all_vcfs.other.subscribe { meta, vcf ->
            if (!excluded_variantcallers.contains(meta.variantcaller)) {
                error("Variant caller '${meta.variantcaller}' is not in the small_variantcallers list and will be excluded from normalization/filtering/consensus. If this is a new SNV caller, please add it to the list in subworkflows/local/post_variantcalling/main.nf:78-80")
            }
        }

        // Needs to be reassigned to enable pass through reassignment below
        // Due to strelka having multiple outputs, we are adding the file name (vcf.gz) for both here to make sure the right files are joined below
        small_variant_vcfs = all_vcfs.small.map{ meta, vcfs_ -> [meta + [filename: vcfs_.name], vcfs_]}
        small_variant_tbis = all_tbis.small.map{ meta, tbis_ -> [meta + [filename: tbis_.baseName], tbis_]}

        // 1. Filter by PASS and custom fields
        // 2. Normalize
        // 3. Aggregate variants (Union, intersection, or n-1)
        if(filter_vcfs) {

            // Join VCFs with their corresponding TBIs before filtering
            FILTER_VCFS( small_variant_vcfs.join(small_variant_tbis, failOnDuplicate: true, failOnMismatch: true), [], [], [])

            small_variant_vcfs = FILTER_VCFS.out.vcf
            small_variant_tbis = FILTER_VCFS.out.tbi
            versions = versions.mix(FILTER_VCFS.out.versions)
        }

        if (normalize_vcfs) {

            NORMALIZE_VCFS(small_variant_vcfs, fasta)

            small_variant_vcfs = NORMALIZE_VCFS.out.vcfs // [meta, vcf]
            small_variant_tbis = NORMALIZE_VCFS.out.tbis // [meta, tbi]
            versions = versions.mix(NORMALIZE_VCFS.out.versions)
        }

        if (normalize_vcfs && snv_consensus_calling){

            CONSENSUS(small_variant_vcfs.join(small_variant_tbis, failOnDuplicate: true, failOnMismatch: true))

            small_variant_vcfs = CONSENSUS.out.vcfs.map { meta, vcfs_ ->
                                        meta.variantcaller = 'consensus'
                                        [meta, vcfs_]
                                    } // [meta, vcfs]

            small_variant_tbis = CONSENSUS.out.tbis.map { meta, tbis_ ->
                                        meta.variantcaller = 'consensus'
                                        [meta, tbis_]
                                    } // [meta, tbis]

            versions = versions.mix(CONSENSUS.out.versions)
        }

        vcfs = small_variant_vcfs.mix(all_vcfs.other)
        tbis = small_variant_tbis.mix(all_tbis.other)

        if (concatenate_vcfs) {
            CONCATENATE_GERMLINE_VCFS(germline_vcfs)

            vcfs = vcfs.mix(CONCATENATE_GERMLINE_VCFS.out.vcfs)
            tbis = tbis.mix(CONCATENATE_GERMLINE_VCFS.out.tbis)

            versions = versions.mix(CONCATENATE_GERMLINE_VCFS.out.versions)
        }


    } else {
        // No post-processing requested, pass through original VCFs
        vcfs = vcfs.mix(germline_vcfs,tumor_only_vcfs, somatic_vcfs)
    }

    emit:
    vcfs     // post processed vcfs [meta, vcf]
    tbis     // post processed tbis [meta, tbi]
    versions // channel: [ versions.yml ]
}
