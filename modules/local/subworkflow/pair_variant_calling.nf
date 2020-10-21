/*
================================================================================
                             SOMATIC VARIANT CALLING
================================================================================
*/

params.strelka_options                = [:]

include { STRELKA_SOMATIC as STRELKA }             from '../../nf-core/software/strelka/somatic'     addParams(options: params.strelka_options)

workflow PAIR_VARIANT_CALLING {
    take:
        bam        // channel: [mandatory] bam
        dbsnp      // channel: [mandatory] dbsnp
        dbsnp_tbi  // channel: [mandatory] dbsnp_tbi
        dict       // channel: [mandatory] dict
        fai        // channel: [mandatory] fai
        fasta      // channel: [mandatory] fasta
        intervals  // channel: [mandatory] intervals
        target_bed // channel: [optional]  target_bed
        tools      //   list:  [mandatory] list of tools

    main:

    bam.map{ meta, bam, bai ->
        patient = meta.patient
        sample  = meta.sample
        gender  = meta.gender
        status  = meta.status
        [patient, sample, gender, status, bam, bai]
    }.branch{
        normal: it[3] == 0
        tumor:  it[3] == 1
    }.set{ bam_to_cross }

    bam_pair = bam_to_cross.normal.cross(bam_to_cross.tumor).map { normal, tumor ->
        def meta = [:]
        meta.patient = normal[0]
        meta.normal  = normal[1]
        meta.tumor   = tumor[1]
        meta.gender  = normal[2]
        meta.id      = "${meta.tumor}_vs_${meta.normal}"

        [meta, normal[4], normal[5], tumor[4], tumor[5]]
    }

    strelka_vcf          = Channel.empty()

    if ('strelka' in tools) {
        STRELKA(
            bam_pair,
            fasta,
            fai,
            target_bed)

        strelka_indels_vcf = STRELKA.out.indels_vcf
        strelka_snvs_vcf   = STRELKA.out.snvs_vcf
        strelka_vcf = strelka_indels_vcf.mix(strelka_snvs_vcf)
    }

    emit:
        strelka_vcf          = strelka_vcf
}
