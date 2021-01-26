/*
================================================================================
                                     MAPPING
================================================================================
*/

params.bwamem1_mem_options    = [:]
params.bwamem2_mem_options    = [:]
params.merge_bam_options      = [:]
params.qualimap_bamqc_options = [:]
params.samtools_index_options = [:]
params.samtools_stats_options = [:]

include { BWA_MEM as BWAMEM1_MEM } from '../process/bwa_mem'                    addParams(options: params.bwamem1_mem_options)
include { BWAMEM2_MEM }            from '../process/bwamem2_mem'                addParams(options: params.bwamem2_mem_options)
include { MERGE_BAM }              from '../process/merge_bam'                  addParams(options: params.merge_bam_options)
include { QUALIMAP_BAMQC }         from '../../nf-core/software/qualimap_bamqc' addParams(options: params.qualimap_bamqc_options)
include { SAMTOOLS_INDEX }         from '../../nf-core/software/samtools/index' addParams(options: params.samtools_index_options)
include { SAMTOOLS_STATS }         from '../../nf-core/software/samtools/stats' addParams(options: params.samtools_stats_options)

workflow MAPPING {
    take:
        skip_bamqc      // boolean: true/false
        skip_samtools   // boolean: true/false
        bwa             // channel: [mandatory] bwa
        fai             // channel: [mandatory] fai
        fasta           // channel: [mandatory] fasta
        reads_input     // channel: [mandatory] reads_input
        save_bam_mapped // boolean: true/false
        step            //   value: [mandatory] starting step
        target_bed      // channel: [optional]  target_bed

    main:

    bam_mapped_index = Channel.empty()
    bam_reports      = Channel.empty()

    if (step == "mapping") {
        bam_bwamem1 = Channel.empty()
        bam_bwamem2 = Channel.empty()

        if (params.aligner == "bwa-mem") {
            BWAMEM1_MEM(reads_input, bwa, fasta, fai)
            bam_bwamem1 = BWAMEM1_MEM.out.bam
        } else {
            BWAMEM2_MEM(reads_input, bwa, fasta, fai)
            bam_bwamem2 = BWAMEM2_MEM.out
        }

        bam_bwa = bam_bwamem1.mix(bam_bwamem2)

        bam_bwa.map{ meta, bam ->
            patient = meta.patient
            sample  = meta.sample
            gender  = meta.gender
            status  = meta.status
            [patient, sample, gender, status, bam]
        }.groupTuple(by: [0,1])
            .branch{
                single:   it[4].size() == 1
                multiple: it[4].size() > 1
            }.set{ bam_bwa_to_sort }

        bam_bwa_single = bam_bwa_to_sort.single.map {
            patient, sample, gender, status, bam ->

            def meta = [:]
            meta.patient = patient
            meta.sample = sample
            meta.gender = gender[0]
            meta.status = status[0]
            meta.id = sample

            [meta, bam[0]]
        }

        bam_bwa_multiple = bam_bwa_to_sort.multiple.map {
            patient, sample, gender, status, bam ->

            def meta = [:]
            meta.patient = patient
            meta.sample = sample
            meta.gender = gender[0]
            meta.status = status[0]
            meta.id = sample

            [meta, bam]
        }

        // STEP 1.5: MERGING AND INDEXING BAM FROM MULTIPLE LANES 
        
        MERGE_BAM(bam_bwa_multiple)
        bam_mapped       = bam_bwa_single.mix(MERGE_BAM.out.bam)
        bam_mapped_index = SAMTOOLS_INDEX(bam_mapped)

        qualimap_bamqc = Channel.empty()
        samtools_stats = Channel.empty()

        if (!skip_bamqc) {
            QUALIMAP_BAMQC(bam_mapped, target_bed)
            qualimap_bamqc = QUALIMAP_BAMQC.out
        }

        if (!skip_samtools) {
            SAMTOOLS_STATS(bam_mapped)
            samtools_stats = SAMTOOLS_STATS.out
        }

        bam_reports = samtools_stats.mix(qualimap_bamqc)

        if (save_bam_mapped) {
            tsv_bam_mapped = bam_mapped.map { meta, bam -> [meta] }
            // Creating TSV files to restart from this step
            tsv_bam_mapped.collectFile(storeDir: "${params.outdir}/preprocessing/tsv") { meta ->
                patient = meta.patient[0]
                sample  = meta.sample[0]
                gender  = meta.gender[0]
                status  = meta.status[0]
                bam   = "${params.outdir}/preprocessing/${sample}/mapped/${sample}.bam"
                bai   = "${params.outdir}/preprocessing/${sample}/mapped/${sample}.bam.bai"
                ["mapped_${sample}.tsv", "${patient}\t${gender}\t${status}\t${sample}\t${bam}\t${bai}\n"]
            }

            tsv_bam_mapped.map { meta ->
                patient = meta.patient[0]
                sample  = meta.sample[0]
                gender  = meta.gender[0]
                status  = meta.status[0]
                bam   = "${params.outdir}/preprocessing/${sample}/mapped/${sample}.bam"
                bai   = "${params.outdir}/preprocessing/${sample}/mapped/${sample}.bam.bai"
                "${patient}\t${gender}\t${status}\t${sample}\t${bam}\t${bai}\n"
            }.collectFile(name: "mapped.tsv", sort: true, storeDir: "${params.outdir}/preprocessing/tsv")
        }
    }

    emit:
        bam = bam_mapped_index
        qc  = bam_reports
}
