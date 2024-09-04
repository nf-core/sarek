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

    ASCAT(cram_pair, allele_files, loci_files, intervals_bed, fasta, gc_file, rt_file)

    ch_versions = ch_versions.mix(ASCAT.out.versions)

    if (asmultipcf) {
        // Group ASCAT outputs by patient
        tumor_logr_by_patient = ASCAT.out.logrs.map { meta, file -> [meta.patient, meta, file] }
            .groupTuple(by: 0)
            .map { patient, metas, files -> [metas[0] + [id:patient], files] }

        tumor_baf_by_patient = ASCAT.out.bafs.map { meta, file -> [meta.patient, meta, file] }
            .groupTuple(by: 0)
            .map { patient, metas, files -> [metas[0] + [id:patient], files] }

        // Assuming normal samples are the same for all tumors of a patient
        normal_logr_by_patient = ASCAT.out.logrs.map { meta, file -> [meta.patient, meta, file] }
            .groupTuple(by: 0)
            .map { patient, metas, files -> [metas[0] + [id:patient], files[0]] }

        normal_baf_by_patient = ASCAT.out.bafs.map { meta, file -> [meta.patient, meta, file] }
            .groupTuple(by: 0)
            .map { patient, metas, files -> [metas[0] + [id:patient], files[0]] }

        // Combine all inputs for ASMULTIPCF
        asmultipcf_input = tumor_logr_by_patient.join(tumor_baf_by_patient)
            .join(normal_logr_by_patient)
            .join(normal_baf_by_patient)

        // Run ASMULTIPCF
        ASMULTIPCF(asmultipcf_input)

        ch_versions = ch_versions.mix(ASMULTIPCF.out.versions)
    }

    emit:
    versions = ch_versions
}
