/*
================================================================================
                                MAPPING
================================================================================
*/

include { BWA_MEM as BWAMEM1_MEM }                  from '../process/bwa_mem'
include { BWAMEM2_MEM }                             from '../process/bwamem2_mem'
include { MERGE_BAM as MERGE_BAM_MAPPED }           from '../process/merge_bam'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_MAPPED } from '../../nf-core/software/samtools/index'

workflow MAPPING {
    take:
        reads_input
        bwa
        fasta
        fai

    main:

    bam_bwamem1 = Channel.empty()
    bam_bwamem2 = Channel.empty()

    if (params.aligner == "bwa-mem") {
        BWAMEM1_MEM(reads_input, bwa, fasta, fai, params.modules['bwa_mem'])
        bam_bwamem1 = BWAMEM1_MEM.out.bam
    } else {
        BWAMEM2_MEM(reads_input, bwa, fasta, fai, params.modules['bwamem2_mem'])
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
    
    MERGE_BAM_MAPPED(bam_bwa_multiple)
    bam_mapped = bam_bwa_single.mix(MERGE_BAM_MAPPED.out.bam)
    bam_mapped = SAMTOOLS_INDEX_MAPPED(bam_mapped, params.modules['samtools_index_mapped'])

    if (params.save_bam_mapped) {
        tsv_bam_mapped = bam_mapped.map { meta, bam, bai -> [meta] }
        // Creating TSV files to restart from this step
        tsv_bam_mapped.collectFile(storeDir: "${params.outdir}/Preprocessing/TSV") { meta ->
            patient = meta.patient[0]
            sample  = meta.sample[0]
            gender  = meta.gender[0]
            status  = meta.status[0]
            bam   = "${params.outdir}/Preprocessing/${sample}/Mapped/${sample}.bam"
            bai   = "${params.outdir}/Preprocessing/${sample}/Mapped/${sample}.bam.bai"
            ["mapped_${sample}.tsv", "${patient}\t${gender}\t${status}\t${sample}\t${bam}\t${bai}\n"]
        }

        tsv_bam_mapped.map { meta ->
            patient = meta.patient[0]
            sample  = meta.sample[0]
            gender  = meta.gender[0]
            status  = meta.status[0]
            bam   = "${params.outdir}/Preprocessing/${sample}/Mapped/${sample}.bam"
            bai   = "${params.outdir}/Preprocessing/${sample}/Mapped/${sample}.bam.bai"
            "${patient}\t${gender}\t${status}\t${sample}\t${bam}\t${bai}\n"
        }.collectFile(name: "mapped.tsv", sort: true, storeDir: "${params.outdir}/Preprocessing/TSV")
    }

    emit:
        bam_mapped = bam_mapped
}
