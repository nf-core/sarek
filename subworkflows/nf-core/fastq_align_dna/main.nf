#!/usr/bin/env nextflow

//
// FASTQ_ALIGN_DNA: Align fastq files to a reference genome
//


include { BOWTIE2_ALIGN                     } from "../../../modules/nf-core/bowtie2/align/main"
include { BWA_MEM as BWAMEM1_MEM            } from '../../../modules/nf-core/bwa/mem/main'
include { BWAMEM2_MEM as BWAMEM2_MEM        } from '../../../modules/nf-core/bwamem2/mem/main'
include { DRAGMAP_ALIGN                     } from "../../../modules/nf-core/dragmap/align/main"
include { SNAPALIGNER_ALIGN as SNAP_ALIGN   } from '../../../modules/nf-core/snapaligner/align/main'



workflow FASTQ_ALIGN_DNA {
    take:
        ch_reads            // channel: [mandatory] meta, reads
        ch_aligner_index    // channel: [mandatory] aligner index
        aligner             // string:  [mandatory] aligner [bowtie2, bwamem, bwamem2, dragmap, snap]
        sort                // boolean: [mandatory] true -> sort, false -> don't sort

    main:

        ch_bai      = Channel.empty()
        ch_bam      = Channel.empty()
        ch_reports  = Channel.empty()
        ch_versions = Channel.empty()

        // Align fastq files to reference genome and (optionally) sort
        switch (aligner) {
            case 'bowtie2':
                BOWTIE2_ALIGN(ch_reads, ch_aligner_index, false, sort) // if aligner is bowtie2
                ch_bam = ch_bam.mix(BOWTIE2_ALIGN.out.bam)
                ch_versions = ch_versions.mix(BOWTIE2_ALIGN.out.versions)
                break
            case 'bwamem':
                BWAMEM1_MEM  (ch_reads, ch_aligner_index, sort)        // If aligner is bwa-mem
                ch_bam = ch_bam.mix(BWAMEM1_MEM.out.bam)
                ch_versions = ch_versions.mix(BWAMEM1_MEM.out.versions)
                break
            case 'bwamem2':
                BWAMEM2_MEM  (ch_reads, ch_aligner_index, sort)        // If aligner is bwa-mem2
                ch_bam = ch_bam.mix(BWAMEM2_MEM.out.bam)
                ch_versions = ch_versions.mix(BWAMEM2_MEM.out.versions)
                break
            case 'dragmap':
                DRAGMAP_ALIGN(ch_reads, ch_aligner_index, sort)        // If aligner is dragmap
                ch_bam = ch_bam.mix(DRAGMAP_ALIGN.out.bam)
                ch_reports = ch_reports.mix(DRAGMAP_ALIGN.out.log)
                ch_versions = ch_versions.mix(DRAGMAP_ALIGN.out.versions)
                break
            case 'snap':
                SNAP_ALIGN   (ch_reads, ch_aligner_index)              // If aligner is snap
                ch_bam = ch_bam.mix(SNAP_ALIGN.out.bam)
                ch_bai.mix(SNAP_ALIGN.out.bai)
                ch_versions = ch_versions.mix(SNAP_ALIGN.out.versions)
                break
            default:
                exit 1, "Unknown aligner: ${aligner}"
        }

    emit:
        bam      = ch_bam         // channel: [ [meta], bam  ]
        bai      = ch_bai         // channel: [ [meta], bai  ]
        reports  = ch_reports     // channel: [ [meta], log  ]
        versions = ch_versions    // channel: [ versions.yml ]
}
