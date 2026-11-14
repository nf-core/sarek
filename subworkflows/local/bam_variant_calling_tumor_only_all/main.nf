//
// TUMOR ONLY VARIANT CALLING
// Should be only run on patients without normal sample
//

include { BAM_VARIANT_CALLING_CNVKIT                             } from '../bam_variant_calling_cnvkit'
include { BAM_VARIANT_CALLING_FREEBAYES                          } from '../bam_variant_calling_freebayes'
include { BAM_VARIANT_CALLING_MPILEUP                            } from '../bam_variant_calling_mpileup'
include { BAM_VARIANT_CALLING_SINGLE_TIDDIT                      } from '../bam_variant_calling_single_tiddit'
include { BAM_VARIANT_CALLING_VARNET                             } from '../bam_variant_calling_varnet'
include { BAM_VARIANT_CALLING_TUMOR_ONLY_CONTROLFREEC            } from '../bam_variant_calling_tumor_only_controlfreec'
include { BAM_VARIANT_CALLING_TUMOR_ONLY_MANTA                   } from '../bam_variant_calling_tumor_only_manta'
include { BAM_VARIANT_CALLING_MUTECT2                            } from '../bam_variant_calling_mutect2'
include { BAM_VARIANT_CALLING_TUMOR_ONLY_LOFREQ                  } from '../bam_variant_calling_tumor_only_lofreq'
include { BAM_VARIANT_CALLING_PARABRICKS_MUTECTCALLER            } from '../bam_variant_calling_parabricks_mutectcaller'
include { BAM_VARIANT_CALLING_TUMOR_ONLY_TNSCOPE                 } from '../bam_variant_calling_tumor_only_tnscope'
include { MSISENSOR2_MSI                                         } from '../../../modules/nf-core/msisensor2/msi'
include { VARNET_GETREF                                          } from '../../../modules/nf-core/varnet/getref'

workflow BAM_VARIANT_CALLING_TUMOR_ONLY_ALL {
    take:
    tools                         // Mandatory, list of tools to apply
    bam                           // channel: [mandatory] bam
    cram                          // channel: [mandatory] cram
    bwa                           // channel: [optional] bwa
    cf_chrom_len                  // channel: [optional] controlfreec length file
    chr_files
    cnvkit_reference
    dbsnp                         // channel: [mandatory] dbsnp
    dbsnp_tbi                     // channel: [mandatory] dbsnp_tbi
    dict                          // channel: [mandatory] dict
    fasta                         // channel: [mandatory] fasta
    fasta_fai                     // channel: [mandatory] fasta_fai
    force_call_alleles            // channel: [optional]  force_call_alleles
    force_call_alleles_tbi        // channel: [optional]  force_call_alleles_tbi
    germline_resource             // channel: [optional]  germline_resource
    germline_resource_tbi         // channel: [optional]  germline_resource_tbi
    contamination_resource        // channel: [optional]  contamination_resource
    contamination_resource_tbi    // channel: [optional]  contamination_resource_tbi
    intervals                     // channel: [mandatory] [ intervals, num_intervals ] or [ [], 0 ] if no intervals
    intervals_bed_combined        // channel: [mandatory] intervals/target regions in one file unzipped
    intervals_bed_gz_tbi_combined // channel: [mandatory] intervals/target regions in one file zipped
    mappability
    msisensor2_models
    panel_of_normals              // channel: [optional]  panel_of_normals
    panel_of_normals_tbi          // channel: [optional]  panel_of_normals_tbi
    varnet_resources              // channel: [optional]  folder with the VarNet germline resources, downloaded if not set
    joint_mutect2                 // boolean: [mandatory] [default: false] run mutect2 in joint mode
    wes                           // boolean: [mandatory] [default: false] whether targeted data is processed
    contamination_resource_provided // boolean: [mandatory] whether a contamination resource was provided

    main:
    // Channels are often remapped to match module/subworkflow

    //TODO: Temporary until the if's can be removed and printing to terminal is prevented with "when" in the modules.config
    out_msisensor2 = channel.empty()
    vcf_freebayes  = channel.empty()
    vcf_lofreq     = channel.empty()
    vcf_manta      = channel.empty()
    vcf_mpileup    = channel.empty()
    vcf_mutect2    = channel.empty()
    vcf_tiddit     = channel.empty()
    vcf_tnscope    = channel.empty()

    // Initialize empty TBI channels
    tbi_freebayes  = channel.empty()
    tbi_lofreq     = channel.empty()
    tbi_manta      = channel.empty()
    tbi_mpileup    = channel.empty()
    tbi_mutect2    = channel.empty()
    tbi_tiddit     = channel.empty()
    tbi_tnscope    = channel.empty()

    // MPILEUP
    if ('mpileup' in tools || 'controlfreec' in tools) {
        BAM_VARIANT_CALLING_MPILEUP(
            cram,
            dict,
            fasta,
            fasta_fai,
            intervals,
            tools,
        )
        vcf_mpileup = BAM_VARIANT_CALLING_MPILEUP.out.vcf
        tbi_mpileup = BAM_VARIANT_CALLING_MPILEUP.out.tbi
    }

    // CONTROLFREEC (depends on MPILEUP)
    if ('controlfreec' in tools) {
        BAM_VARIANT_CALLING_TUMOR_ONLY_CONTROLFREEC(
            BAM_VARIANT_CALLING_MPILEUP.out.mpileup.map { meta, pileup_tumor -> [meta, [], pileup_tumor, [], [], [], []] },
            fasta.map { _meta, fasta_ -> [fasta_] },
            cf_chrom_len ?: fasta_fai.map { _meta, fasta_fai_ -> [fasta_fai_] },
            dbsnp,
            dbsnp_tbi,
            chr_files,
            mappability,
            wes ? intervals_bed_combined : [],
        )

    }

    // CNVKIT
    if ('cnvkit' in tools) {
        BAM_VARIANT_CALLING_CNVKIT(
            bam.map { meta_, bam_, _bai -> [meta_, bam_, []] },
            fasta,
            fasta_fai,
            [[id: "null"], []],
            cnvkit_reference.map { reference -> [[id: reference[0].baseName], reference] },
        )
    }

    // FREEBAYES
    if ('freebayes' in tools) {
        BAM_VARIANT_CALLING_FREEBAYES(
            cram.map { meta_, cram_, crai_ -> [meta_, cram_, crai_, [], []] },
            dict,
            fasta,
            fasta_fai,
            intervals,
        )

        vcf_freebayes = BAM_VARIANT_CALLING_FREEBAYES.out.vcf
        tbi_freebayes = BAM_VARIANT_CALLING_FREEBAYES.out.tbi
    }

    // MSISENSOR
    if ('msisensor2' in tools) {

        MSISENSOR2_MSI(bam, msisensor2_models)

        out_msisensor2 = out_msisensor2.mix(MSISENSOR2_MSI.out.distribution)
        out_msisensor2 = out_msisensor2.mix(MSISENSOR2_MSI.out.somatic)
    }

    // Mutect2 and Parabricks MutectCaller take [ meta, normal_cram, normal_crai, tumor_cram, tumor_crai ], with no normal here
    cram_tumor_only = cram.map { meta, tumor_cram, tumor_crai -> [meta, [], [], tumor_cram, tumor_crai] }

    // MUTECT2
    if ('mutect2' in tools) {
        BAM_VARIANT_CALLING_MUTECT2(
            cram_tumor_only,
            fasta,
            fasta_fai,
            dict,
            force_call_alleles,
            force_call_alleles_tbi,
            germline_resource,
            germline_resource_tbi,
            contamination_resource,
            contamination_resource_tbi,
            panel_of_normals,
            panel_of_normals_tbi,
            intervals,
            joint_mutect2,
            contamination_resource_provided,
        )

        vcf_mutect2 = BAM_VARIANT_CALLING_MUTECT2.out.vcf
        tbi_mutect2 = BAM_VARIANT_CALLING_MUTECT2.out.tbi
    }

    // PARABRICKS MUTECTCALLER
    cram_parabricks_mutectcaller = cram_tumor_only.filter { _meta, _normal_cram, _normal_crai, _tumor_cram, _tumor_crai ->
        'parabricks_mutectcaller' in tools
    }

    BAM_VARIANT_CALLING_PARABRICKS_MUTECTCALLER(
        cram_parabricks_mutectcaller,
        fasta,
        fasta_fai,
        dict,
        force_call_alleles,
        force_call_alleles_tbi,
        germline_resource,
        germline_resource_tbi,
        contamination_resource,
        contamination_resource_tbi,
        panel_of_normals,
        panel_of_normals_tbi,
        intervals_bed_combined,
        intervals,
        contamination_resource_provided,
    )

    vcf_parabricks_mutectcaller = BAM_VARIANT_CALLING_PARABRICKS_MUTECTCALLER.out.vcf
    tbi_parabricks_mutectcaller = BAM_VARIANT_CALLING_PARABRICKS_MUTECTCALLER.out.tbi

    //LOFREQ
    if ('lofreq' in tools) {
        BAM_VARIANT_CALLING_TUMOR_ONLY_LOFREQ(
            cram,
            fasta,
            fasta_fai,
            intervals,
            dict,
        )
        vcf_lofreq = BAM_VARIANT_CALLING_TUMOR_ONLY_LOFREQ.out.vcf
        tbi_lofreq = BAM_VARIANT_CALLING_TUMOR_ONLY_LOFREQ.out.tbi
    }

    // MANTA
    if ('manta' in tools) {
        BAM_VARIANT_CALLING_TUMOR_ONLY_MANTA(
            cram,
            fasta,
            fasta_fai,
            intervals_bed_gz_tbi_combined,
        )

        vcf_manta = BAM_VARIANT_CALLING_TUMOR_ONLY_MANTA.out.tumor_sv_vcf
        tbi_manta = BAM_VARIANT_CALLING_TUMOR_ONLY_MANTA.out.tumor_sv_vcf_tbi
    }

    // TIDDIT
    if ('tiddit' in tools) {
        BAM_VARIANT_CALLING_SINGLE_TIDDIT(
            cram,
            fasta,
            fasta_fai,
            bwa,
        )

        vcf_tiddit = BAM_VARIANT_CALLING_SINGLE_TIDDIT.out.vcf
        tbi_tiddit = BAM_VARIANT_CALLING_SINGLE_TIDDIT.out.tbi
    }

    // TNSCOPE
    if ('sentieon_tnscope' in tools) {
        BAM_VARIANT_CALLING_TUMOR_ONLY_TNSCOPE(
            cram,
            fasta,
            fasta_fai,
            dict,
            germline_resource,
            germline_resource_tbi,
            panel_of_normals,
            panel_of_normals_tbi,
            intervals,
        )

        vcf_tnscope = BAM_VARIANT_CALLING_TUMOR_ONLY_TNSCOPE.out.vcf
        tbi_tnscope = BAM_VARIANT_CALLING_TUMOR_ONLY_TNSCOPE.out.tbi
    }

    // VARNET
    bam_varnet = bam
        .filter { _meta, _bam, _bai -> 'varnet' in tools }
        .map { meta_, bam_, bai -> [meta_, [], [], bam_, bai] }

    // download the germline resources once if --varnet_resources is not set and VarNet runs on tumor-only samples
    VARNET_GETREF(
        varnet_resources
            .filter { dir -> !dir }
            .combine(bam_varnet.first())
            .map { _resources_and_bam -> [id: 'varnet_resources'] }
    )

    varnet_dbsnp = varnet_resources
        .filter { dir -> dir }
        .map { dir -> [[id: 'varnet_resources'], file("${dir}/GrCh38-dbSNP-build-156-common-variants-0-indexed-CHROMPOS.npz", checkIfExists: true)] }
        .mix(VARNET_GETREF.out.dbsnp)
        .first()
    varnet_gnomad = varnet_resources
        .filter { dir -> dir }
        .map { dir -> [[id: 'varnet_resources'], file("${dir}/combined_gnomAD_v4.1_genomes_exomes_dict.hg38.chunks.npz", checkIfExists: true)] }
        .mix(VARNET_GETREF.out.gnomad)
        .first()
    varnet_pon = varnet_resources
        .filter { dir -> dir }
        .map { dir -> [[id: 'varnet_resources'], file("${dir}/1000g_pon.hg38.vcf.gz", checkIfExists: true)] }
        .mix(VARNET_GETREF.out.pon)
        .first()

    BAM_VARIANT_CALLING_VARNET(
        bam_varnet,
        fasta,
        fasta_fai,
        intervals_bed_combined,
        varnet_dbsnp,
        varnet_gnomad,
        varnet_pon,
    )

    vcf_varnet = BAM_VARIANT_CALLING_VARNET.out.vcf
    tbi_varnet = BAM_VARIANT_CALLING_VARNET.out.tbi

    vcf_all = channel.empty()
        .mix(
            vcf_freebayes,
            vcf_lofreq,
            vcf_manta,
            vcf_mutect2,
            vcf_mpileup,
            vcf_parabricks_mutectcaller,
            vcf_tiddit,
            vcf_tnscope,
            vcf_varnet,
        )

    tbi_all = channel.empty()
        .mix(
            tbi_freebayes,
            tbi_lofreq,
            tbi_manta,
            tbi_mutect2,
            tbi_mpileup,
            tbi_parabricks_mutectcaller,
            tbi_tiddit,
            tbi_tnscope,
            tbi_varnet,
        )

    emit:
    out_msisensor2
    vcf_all
    tbi_all
    vcf_freebayes
    vcf_lofreq
    vcf_manta
    vcf_mpileup
    vcf_mutect2
    vcf_parabricks_mutectcaller
    vcf_tiddit
    vcf_tnscope
    vcf_varnet
    tbi_freebayes
    tbi_lofreq
    tbi_manta
    tbi_mpileup
    tbi_mutect2
    tbi_parabricks_mutectcaller
    tbi_tiddit
    tbi_tnscope
    tbi_varnet
}
