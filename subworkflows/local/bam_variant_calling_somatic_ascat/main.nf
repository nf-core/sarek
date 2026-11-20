//
// ASCAT variant calling
//
// For all modules here:
// A when clause condition is defined in the conf/modules.config to determine if the module should be run

include { ASCAT } from '../../../modules/nf-core/ascat/main'

workflow BAM_VARIANT_CALLING_SOMATIC_ASCAT {

    take:
    cram_pair                // channel: [mandatory] [meta, normal_cram, normal_crai, tumor_cram, tumor_crai]
    allele_files             // channel: [mandatory] zip
    loci_files               // channel: [mandatory] zip
    intervals_bed            // channel: [optional]  bed for WES
    fasta                    // channel: [optional]  fasta needed for cram
    gc_file                  // channel: [optional]  txt for LogRCorrection
    rt_file                  // channel: [optional]  txt for LogRCorrection
    ascat_genome             // value:   [optional]  genome build (hg19 or hg38), read by ext.args through meta

    main:

    // ext.args cannot read a workflow value, so the genome build travels in meta
    cram_pair_with_genome = cram_pair.map { meta, normal_cram, normal_crai, tumor_cram, tumor_crai ->
        [meta + [ascat_genome: ascat_genome], normal_cram, normal_crai, tumor_cram, tumor_crai]
    }

    ASCAT(cram_pair_with_genome, allele_files, loci_files, intervals_bed, fasta, gc_file, rt_file)
}
