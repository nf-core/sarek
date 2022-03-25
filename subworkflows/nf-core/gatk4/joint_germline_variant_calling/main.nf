//
// Run GATK haplotypecaller for all input samples, merge them with genomicsdbimport, perform joint genotyping with genotypeGVCFS and recalibrate with variantrecalibrator & applyvqsr.
//

include { GATK4_HAPLOTYPECALLER     as HAPLOTYPECALLER     } from '../../../../modules/nf-core/modules/gatk4/haplotypecaller/main'
include { GATK4_GENOMICSDBIMPORT    as GENOMICSDBIMPORT    } from '../../../../modules/nf-core/modules/gatk4/genomicsdbimport/main'
include { GATK4_GENOTYPEGVCFS       as GENOTYPEGVCFS       } from '../../../../modules/nf-core/modules/gatk4/genotypegvcfs/main'
include { GATK4_VARIANTRECALIBRATOR as VARIANTRECALIBRATOR } from '../../../../modules/nf-core/modules/gatk4/variantrecalibrator/main'
include { GATK4_APPLYVQSR           as APPLYVQSR           } from '../../../../modules/nf-core/modules/gatk4/applyvqsr/main'

workflow GATK_JOINT_GERMLINE_VARIANT_CALLING {
    take:
    input                     // channel: [ val(meta), [ input ], [ input_index ], [] ]
    fasta                     // channel: /path/to/reference/fasta
    fai                       // channel: /path/to/reference/fasta/index
    dict                      // channel: /path/to/reference/fasta/dictionary
    sites                     // channel: /path/to/known/sites/file
    sites_index               // channel: /path/to/known/sites/index
    allelespecific            // channel: true/false run allelespecific mode of vqsr modules
    resources                 // channel: [[resource, vcfs, forvariantrecal], [resource, tbis, forvariantrecal], [resource, labels, forvariantrecal]]
    annotation                // channel: [annotations, to, use, for, variantrecal, filtering]
    mode                      // channel: which mode to run variantrecal: SNP/INDEL/BOTH
    create_rscript            // channel: true/false whether to generate rscript plots in variantrecal
    truthsensitivity          // channel: 0-100.0 truthsensitivity cutoff for applyvqsr

    main:
    ch_versions       = Channel.empty()

    //
    //Convert all sample vcfs into a genomicsdb workspace using genomicsdbimport.
    //
    gendb_input       = .combine([interval_file]).combine(['']).combine([dict])

    GENOMICSDBIMPORT ( gendb_input, false, false, false )
    ch_versions       = ch_versions.mix(GENOMICSDBIMPORT.out.versions)

    //
    //Joint genotyping performed using GenotypeGVCFs
    //
    // ch_genotype_in    = GENOMICSDBIMPORT.out.genomicsdb.collect()
    // //[] is a placeholder for the input where the vcf tbi file would be passed in for non-genomicsdb workspace runs, which is not necessary for this workflow as it uses genomicsdb workspaces.
    // ch_genotype_in.add([])

    // GENOTYPEGVCFS ( ch_genotype_in, fasta, fai, dict, sites, sites_index )
    // ch_versions       = ch_versions.mix(GENOTYPEGVCFS.out.versions)

    // // setting run_vqsr to false skips the VQSR process, for if user does not wish to perform VQSR,
    // // or want to run the hard filtering recommended by gatk best practices for runs with a low number of samples instead.
    // if (run_vqsr) {
    //     //
    //     //Perform first step in VQSR using VariantRecalibrator
    //     //
    //     ch_gvcf       = GENOTYPEGVCFS.out.vcf.collect()
    //     ch_gtbi       = GENOTYPEGVCFS.out.tbi.collect()
    //     ch_vrecal_in  = ch_gvcf.combine(ch_gtbi, by: 0)

    //     VARIANTRECALIBRATOR ( ch_vrecal_in, fasta, fai, dict, allelespecific, resources, annotation, mode, create_rscript )

    //     ch_versions   = ch_versions.mix(VARIANTRECALIBRATOR.out.versions)

    //     //
    //     //Perform second step in VQSR using ApplyVQSR
    //     //
    //     ch_recal      = VARIANTRECALIBRATOR.out.recal.collect()
    //     ch_idx        = VARIANTRECALIBRATOR.out.idx.collect()
    //     ch_tranches   = VARIANTRECALIBRATOR.out.tranches.collect()
    //     ch_vqsr_in    = ch_vrecal_in.combine(ch_recal, by: 0).combine(ch_idx, by: 0).combine(ch_tranches, by: 0)

    //     APPLYVQSR ( ch_vqsr_in, fasta, fai, dict, allelespecific, truthsensitivity, mode )

    //     ch_versions   = ch_versions.mix(APPLYVQSR.out.versions)

    // }

    emit:
    // haplotc_vcf    = run_haplotc ? HAPLOTYPECALLER.out.vcf.collect()          : [] // channel: [ val(meta), [ vcf ] ]
    // haplotc_index  = run_haplotc ? HAPLOTYPECALLER.out.tbi.collect()          : [] // channel: [ val(meta), [ tbi ] ]

    // genomicsdb     =               GENOMICSDBIMPORT.out.genomicsdb.collect()       // channel: [ val(meta), [ genomicsdb ] ]

    // genotype_vcf   =               GENOTYPEGVCFS.out.vcf.collect()                 // channel: [ val(meta), [ vcf ] ]
    // genotype_index =               GENOTYPEGVCFS.out.vcf.collect()                 // channel: [ val(meta), [ tbi ] ]

    // recal_file     = run_vqsr    ? VARIANTRECALIBRATOR.out.recal.collect()    : [] // channel: [ val(meta), [ recal ] ]
    // recal_index    = run_vqsr    ? VARIANTRECALIBRATOR.out.idx.collect()      : [] // channel: [ val(meta), [ idx ] ]
    // recal_tranches = run_vqsr    ? VARIANTRECALIBRATOR.out.tranches.collect() : [] // channel: [ val(meta), [ tranches ] ]

    // vqsr_vcf       = run_vqsr    ? APPLYVQSR.out.vcf.collect()                : [] // channel: [ val(meta), [ vcf ] ]
    // vqsr_index     = run_vqsr    ? APPLYVQSR.out.tbi.collect()                : [] // channel: [ val(meta), [ tbi ] ]

    versions       =               ch_versions                                     // channel: [ versions.yml ]
}
