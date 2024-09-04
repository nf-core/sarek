//
// ASCAT variant calling
//
// For all modules here:
// A when clause condition is defined in the conf/modules.config to determine if the module should be run

include { ASCAT } from '../../../modules/nf-core/ascat/main'
include { ASMULTIPCF } from '../../../modules/nf-core/asmultipcf/main'

workflow BAM_VARIANT_CALLING_SOMATIC_ASCAT {
    take:
    cram_pair                // channel: [mandatory] [meta, normal_cram, normal_crai, tumor_cram, tumor_crai]
    allele_files             // channel: [mandatory] zip
    loci_files               // channel: [mandatory] zip
    intervals_bed            // channel: [optional]  bed for WES
    fasta                    // channel: [optional]  fasta needed for cram
    gc_file                  // channel: [optional]  txt for LogRCorrection
    rt_file                  // channel: [optional]  txt for LogRCorrection
    asmultipcf               // boolean: [mandatory] whether to run ASMULTIPCF

    main:
    ch_versions = Channel.empty()

    // Group input by patient
    cram_pair_by_patient = cram_pair
        .map { meta, normal_cram, normal_crai, tumor_cram, tumor_crai -> 
            [meta.patient, meta, normal_cram, normal_crai, tumor_cram, tumor_crai] 
        }
        .groupTuple()

    // Run ASCAT for all samples
    ASCAT(cram_pair, allele_files, loci_files, intervals_bed, fasta, gc_file, rt_file)

    // Group ASCAT outputs by patient
    ascat_output_by_patient = ASCAT.out.logrs.join(ASCAT.out.bafs)
        .map { meta, logr, baf -> [meta.patient, meta, logr, baf] }
        .groupTuple()

    if (asmultipcf) {
        // Prepare input for ASMULTIPCF
        asmultipcf_input = ascat_output_by_patient
            .map { patient, metas, logrs, bafs ->
                def tumor_logrs = logrs.findAll { it.name.contains('tumor') }
                def tumor_bafs = bafs.findAll { it.name.contains('tumor') }
                def normal_logr = logrs.find { it.name.contains('normal') }
                def normal_baf = bafs.find { it.name.contains('normal') }
                [metas[0] + [id: patient], tumor_logrs, tumor_bafs, normal_logr, normal_baf]
            }

        // Run ASMULTIPCF
        ASMULTIPCF(asmultipcf_input)

        ch_versions = ch_versions.mix(ASMULTIPCF.out.versions)
    }

    ch_versions = ch_versions.mix(ASCAT.out.versions)

    emit:
    ascat_segments = ASCAT.out.segments
    ascat_purityploidy = ASCAT.out.purityploidy
    asmultipcf_segments = asmultipcf ? ASMULTIPCF.out.asmultipcf_segments : Channel.empty()
    asmultipcf_purityploidy = asmultipcf ? ASMULTIPCF.out.asmultipcf_purityploidy : Channel.empty()
    versions = ch_versions
}
