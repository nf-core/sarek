//
// Runs FGBIO tools to remove UMI tags from FASTQ reads
// Convert them to unmapped BAM file, map them to the reference genome,
// use the mapped information to group UMIs and generate consensus reads
//
// For all modules here:
// A when clause condition is defined in the conf/modules.config to determine if the module should be run

include { FGBIO_CALLMOLECULARCONSENSUSREADS as CALLUMICONSENSUS } from '../../../modules/nf-core/fgbio/callmolecularconsensusreads/main.nf'
include { FGBIO_FASTQTOBAM                  as FASTQTOBAM       } from '../../../modules/nf-core/fgbio/fastqtobam/main'
include { FGBIO_GROUPREADSBYUMI             as GROUPREADSBYUMI  } from '../../../modules/nf-core/fgbio/groupreadsbyumi/main'
include { FASTQ_ALIGN as ALIGN_UMI                              } from '../fastq_align/main'
include { SAMTOOLS_MERGE                                        } from '../../../modules/nf-core/samtools/merge/main'
include { SAMTOOLS_BAM2FQ                   as BAM2FASTQ        } from '../../../modules/nf-core/samtools/bam2fq/main.nf'

workflow FASTQ_CREATE_UMI_CONSENSUS_FGBIO {
    take:
    reads                     // channel: [mandatory] [ val(meta), [ reads ] ]
    fasta                     // channel: [mandatory] /path/to/reference/fasta
    fai                       // channel: [optional] /path/to/reference/fasta_fai, needed for Sentieon
    map_index                 // channel: [mandatory] Pre-computed mapping index
    groupreadsbyumi_strategy  // string:  [mandatory] grouping strategy - default: "Adjacency"

    main:
    ch_versions = Channel.empty()

    // params.umi_read_structure is passed out as ext.args
    // FASTQ reads are converted into a tagged unmapped BAM file (uBAM)
    FASTQTOBAM(reads)

    // in order to map uBAM using BWA MEM, we need to convert uBAM to FASTQ
    // TODO check if DRAGMAP works well with BAM inputs
    // but keep the appropriate UMI tags in the FASTQ comment field and produce
    // an interleaved FASQT file (hence, split = false)
    split = false
    BAM2FASTQ(FASTQTOBAM.out.bam, split)

    // appropriately tagged interleaved FASTQ reads are mapped to the reference
    // bams will not be sorted (hence, sort = false)
    sort = false
    ALIGN_UMI(BAM2FASTQ.out.reads, map_index, sort, fasta, fai)

    bams_to_merge = ALIGN_UMI.out.bam
    // id currently includes the lane, so swap to just id=sample and groupKey to avoid blocking
        .map {meta, bam ->
            tuple( groupKey(meta + [id:meta.sample], meta.num_lanes), bam)
            }
        .groupTuple()
        // undo the groupKey, else the meta map is not a normal map.
        .map{meta, bam -> tuple(meta.target, bam)}
        .branch { meta, bam ->
            single: meta.num_lanes <= 1
            return [meta, bam[0]]
            multiple: meta.num_lanes > 1
        }

    // Merge across runs/lanes for the same sample
    SAMTOOLS_MERGE(bams_to_merge.multiple, [[], []], [[], []])

    bams_all = SAMTOOLS_MERGE.out.bam.mix(bams_to_merge.single)

    // appropriately tagged reads are now grouped by UMI information
    GROUPREADSBYUMI(bams_all, groupreadsbyumi_strategy)

    // Using newly created groups
    // To call a consensus across reads in the same group
    // And emit a consensus BAM file
    // TODO: add params for call_min_reads and call_min_baseq
    call_min_reads = 1
    call_min_baseq = 10
    CALLUMICONSENSUS(GROUPREADSBYUMI.out.bam, call_min_reads, call_min_baseq)

    ch_versions = ch_versions.mix(BAM2FASTQ.out.versions)
    ch_versions = ch_versions.mix(ALIGN_UMI.out.versions)
    ch_versions = ch_versions.mix(CALLUMICONSENSUS.out.versions)
    ch_versions = ch_versions.mix(FASTQTOBAM.out.versions)
    ch_versions = ch_versions.mix(GROUPREADSBYUMI.out.versions)
    ch_versions = ch_versions.mix(SAMTOOLS_MERGE.out.versions)

    emit:
    umibam         = FASTQTOBAM.out.bam             // channel: [ val(meta), [ bam ] ]
    groupbam       = GROUPREADSBYUMI.out.bam        // channel: [ val(meta), [ bam ] ]
    consensusbam   = CALLUMICONSENSUS.out.bam       // channel: [ val(meta), [ bam ] ]
    versions       = ch_versions                    // channel: [ versions.yml ]
}
