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
    varlociraptor_events_germline
    varlociraptor_events_somatic
    varlociraptor_events_tumor_only
    varlociraptor_fdr
    varlociraptor_scenario_germline
    varlociraptor_scenario_somatic
    varlociraptor_scenario_tumor_only

    main:
    // Post-processed VCFs are carried as [meta, vcf, tbi] and only split into vcfs/tbis at emit
    vcf_tbi          = channel.empty()
    // Original VCFs when no post-processing is requested (passed through without their tbis)
    passthrough_vcfs = channel.empty()

    //
    // VARLOCIRAPTOR
    //
    if (tools && tools.split(',').contains('varlociraptor')) {
        // GERMLINE
        VCF_VARLOCIRAPTOR_GERMLINE(cram_germline, fasta, fai, varlociraptor_scenario_germline, germline_vcfs, varlociraptor_chunk_size, 'normal', varlociraptor_events_germline, varlociraptor_fdr)

        // SOMATIC
        VCF_VARLOCIRAPTOR_SOMATIC(cram_somatic, fasta, fai, varlociraptor_scenario_somatic, somatic_vcfs, germline_vcfs, varlociraptor_chunk_size, varlociraptor_events_somatic, varlociraptor_fdr)

        // TUMOR ONLY
        VCF_VARLOCIRAPTOR_TUMOR_ONLY(cram_tumor_only, fasta, fai, varlociraptor_scenario_tumor_only, tumor_only_vcfs, varlociraptor_chunk_size, 'tumor', varlociraptor_events_tumor_only, varlociraptor_fdr)

        varlociraptor_vcfs = channel.empty().mix(VCF_VARLOCIRAPTOR_GERMLINE.out.vcf, VCF_VARLOCIRAPTOR_SOMATIC.out.vcf, VCF_VARLOCIRAPTOR_TUMOR_ONLY.out.vcf)
        varlociraptor_tbis = channel.empty().mix(VCF_VARLOCIRAPTOR_GERMLINE.out.tbi, VCF_VARLOCIRAPTOR_SOMATIC.out.tbi, VCF_VARLOCIRAPTOR_TUMOR_ONLY.out.tbi)

        vcf_tbi = varlociraptor_vcfs.join(varlociraptor_tbis, failOnDuplicate: true, failOnMismatch: true)
            .map { meta, vcf, tbi -> [meta + [publish_dir: 'varlociraptor'], vcf, tbi] }

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
        def small_variantcallers = ['bcftools', 'deepsomatic', 'deepvariant', 'freebayes', 'haplotypecaller',
                                    'lofreq', 'muse', 'mutect2', 'parabricks_deepsomatic', 'parabricks_deepvariant', 'parabricks_haplotypecaller',
                                    'parabricks_mutectcaller', 'sentieon_dnascope', 'sentieon_haplotyper', 'sentieon_tnscope', 'strelka', 'varnet' ]

        def excluded_variantcallers = ['manta', 'tiddit', 'samtools']

        // Pair each VCF with its tbi. Somatic Strelka and Manta emit two VCFs per meta, so the VCF file name is part of the join key
        all_vcfs = channel.empty().mix(germline_vcfs, tumor_only_vcfs, somatic_vcfs).map { meta, vcf -> [[meta, vcf.name], vcf] }
        all_tbis = channel.empty().mix(germline_tbis, tumor_only_tbis, somatic_tbis).map { meta, tbi -> [[meta, tbi.baseName], tbi] }

        all_vcf_tbi = all_vcfs.join(all_tbis, failOnDuplicate: true, failOnMismatch: true)
                                .map { key, vcf, tbi -> [key[0], vcf, tbi] }
                                .branch{ meta, _vcf, _tbi ->
                                    small: small_variantcallers.contains(meta.variantcaller)
                                    other: true
                                }

        // Validate that we're not silently excluding unknown variant callers
        all_vcf_tbi.other.subscribe { meta, _vcf, _tbi ->
            if (!excluded_variantcallers.contains(meta.variantcaller)) {
                error("Variant caller '${meta.variantcaller}' is not in the small_variantcallers list and will be excluded from normalization/filtering/consensus. If this is a new SNV caller, please add it to the list in subworkflows/local/post_variantcalling/main.nf")
            }
        }

        // Due to somatic strelka having two VCFs per meta (SNVs and indels), the file name is added to meta to keep them apart in the joins below
        small_variant_vcf_tbi = all_vcf_tbi.small.map{ meta, vcf, tbi -> [meta + [filename: vcf.name], vcf, tbi] }
        def small_variant_publish_dir = null

        // 1. Filter by PASS and custom fields
        // 2. Normalize
        // 3. Aggregate variants (Union, intersection, or n-1)
        if(filter_vcfs) {

            FILTER_VCFS(small_variant_vcf_tbi, [], [], [])

            small_variant_vcf_tbi = FILTER_VCFS.out.vcf.join(FILTER_VCFS.out.index, failOnDuplicate: true, failOnMismatch: true)
            small_variant_publish_dir = 'filtered'
        }

        if (normalize_vcfs) {

            NORMALIZE_VCFS(small_variant_vcf_tbi.map { meta, vcf, _tbi -> [meta, vcf] }, fasta)

            small_variant_vcf_tbi = NORMALIZE_VCFS.out.vcfs.join(NORMALIZE_VCFS.out.tbis, failOnDuplicate: true, failOnMismatch: true)
            small_variant_publish_dir = 'normalized'
        }

        // Individual caller VCFs are kept next to any consensus VCF for downstream annotation
        vcf_tbi = small_variant_vcf_tbi
            .map { meta, vcf, tbi -> [meta + [publish_dir: small_variant_publish_dir], vcf, tbi] }
            .mix(all_vcf_tbi.other)

        if (normalize_vcfs && snv_consensus_calling){

            CONSENSUS(small_variant_vcf_tbi)

            // variantcaller 'consensus': published to variant_calling/consensus/, so no meta.publish_dir needed
            vcf_tbi = vcf_tbi.mix(
                CONSENSUS.out.vcfs.join(CONSENSUS.out.tbis, failOnDuplicate: true, failOnMismatch: true)
                    .map { meta, vcf, tbi -> [meta + [variantcaller: 'consensus'], vcf, tbi] }
            )
        }

        if (concatenate_vcfs) {
            CONCATENATE_GERMLINE_VCFS(germline_vcfs)

            vcf_tbi = vcf_tbi.mix(
                CONCATENATE_GERMLINE_VCFS.out.vcfs.join(CONCATENATE_GERMLINE_VCFS.out.tbis, failOnDuplicate: true, failOnMismatch: true)
                    .map { meta, vcf, tbi -> [meta + [publish_dir: 'concat'], vcf, tbi] }
            )
        }

    } else {
        // No post-processing requested, pass through original VCFs
        passthrough_vcfs = channel.empty().mix(germline_vcfs, tumor_only_vcfs, somatic_vcfs)
    }

    emit:
    // meta.publish_dir: variant_calling/ subdirectory a post-processed VCF is published to by conf/modules
    // (post_variant_calling.config, varlociraptor.config); unset means variant_calling/<variantcaller>/.
    // Only read by CHANNEL_VARIANT_CALLING_CREATE_CSV, to point csv/variantcalled.csv at the published file.
    vcfs = vcf_tbi.map { meta, vcf, _tbi -> [meta, vcf] }.mix(passthrough_vcfs) // [meta, vcf]
    tbis = vcf_tbi.map { meta, _vcf, tbi -> [meta, tbi] }                       // [meta, tbi], empty without post-processing
}
