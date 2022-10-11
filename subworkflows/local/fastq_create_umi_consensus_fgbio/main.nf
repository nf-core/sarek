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
include { FASTQ_ALIGN_BWAMEM_MEM2_DRAGMAP   as ALIGN_UMI        } from '../fastq_align_bwamem_mem2_dragmap/main'
include { SAMBLASTER                                            } from '../../../modules/nf-core/samblaster/main'
include { SAMTOOLS_BAM2FQ                   as BAM2FASTQ        } from '../../../modules/nf-core/samtools/bam2fq/main.nf'

workflow FASTQ_CREATE_UMI_CONSENSUS_FGBIO {
    take:
    reads                     // channel: [mandatory] [ val(meta), [ reads ] ]
    fasta                     // channel: [mandatory] /path/to/reference/fasta
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
    ALIGN_UMI(BAM2FASTQ.out.reads, map_index, sort)

    // samblaster is used in order to tag mates information in the BAM file
    // this is used in order to group reads by UMI
    SAMBLASTER(ALIGN_UMI.out.bam)

    // appropriately tagged reads are now grouped by UMI information
    GROUPREADSBYUMI(SAMBLASTER.out.bam, groupreadsbyumi_strategy)

    // Using newly created groups
    // To call a consensus across reads in the same group
    // And emit a consensus BAM file
    CALLUMICONSENSUS(GROUPREADSBYUMI.out.bam)

    ch_versions = ch_versions.mix(BAM2FASTQ.out.versions)
    ch_versions = ch_versions.mix(ALIGN_UMI.out.versions)
    ch_versions = ch_versions.mix(CALLUMICONSENSUS.out.versions)
    ch_versions = ch_versions.mix(FASTQTOBAM.out.versions)
    ch_versions = ch_versions.mix(GROUPREADSBYUMI.out.versions)
    ch_versions = ch_versions.mix(SAMBLASTER.out.versions)

    emit:
    umibam         = FASTQTOBAM.out.bam             // channel: [ val(meta), [ bam ] ]
    groupbam       = GROUPREADSBYUMI.out.bam        // channel: [ val(meta), [ bam ] ]
    consensusbam   = CALLUMICONSENSUS.out.bam       // channel: [ val(meta), [ bam ] ]
    versions       = ch_versions                    // channel: [ versions.yml ]
}
