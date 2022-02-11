//
// Runs FGBIO tools to remove UMI tags from FASTQ reads
// Convert them to unmapped BAM file, map them to the reference genome,
// use the mapped information to group UMIs and generate consensus reads
//


include { BWAMEM2_MEM                                             } from '../../../modules/nf-core/modules/bwamem2/mem/main'
include { BWA_MEM                           as BWAMEM1_MEM        } from '../../../modules/nf-core/modules/bwa/mem/main'
include { FGBIO_CALLMOLECULARCONSENSUSREADS as CALLUMICONSENSUS   } from '../../../modules/nf-core/modules/fgbio/callmolecularconsensusreads/main.nf'
include { FGBIO_FASTQTOBAM                  as FASTQTOBAM         } from '../../../modules/nf-core/modules/fgbio/fastqtobam/main'
include { FGBIO_GROUPREADSBYUMI             as GROUPREADSBYUMI    } from '../../../modules/nf-core/modules/fgbio/groupreadsbyumi/main'
include { SAMBLASTER                                              } from '../../../modules/nf-core/modules/samblaster/main'
include { SAMTOOLS_BAM2FQ                   as BAM2FASTQ          } from '../../../modules/nf-core/modules/samtools/bam2fq/main.nf'

workflow CREATE_UMI_CONSENSUS {
    take:
    reads                     // channel: [mandatory] [ val(meta), [ reads ] ]
    fasta                     // channel: [mandatory] /path/to/reference/fasta
    bwa                       // channel: [mandatory] Pre-computed BWA index (either bwa-mem or bwa-mem2; MUST be matching to chosen aligner)
    read_structure            // string:  [mandatory] "read_structure"
    groupreadsbyumi_strategy  // string:  [mandatory] grouping strategy - default: "Adjacency"
    aligner                   // string:  [mandatory] "bwa-mem" or "bwa-mem2"

    main:
    ch_versions = Channel.empty()

    // using information in val(read_structure) FASTQ reads are converted into
    // a tagged unmapped BAM file (uBAM)
    FASTQTOBAM ( reads, read_structure )
    ch_versions = ch_versions.mix(FASTQTOBAM.out.versions)

    // in order to map uBAM using BWA MEM, we need to convert uBAM to FASTQ
    // but keep the appropriate UMI tags in the FASTQ comment field and produce
    // an interleaved FASQT file (hence, split = false)
    split = false
    BAM2FASTQ ( FASTQTOBAM.out.umibam, split )
    ch_versions = ch_versions.mix(BAM2FASTQ.out.versions)

    // the user can choose here to use either bwa-mem (default) or bwa-mem2
    aligned_bam = Channel.empty()

    if (aligner == "bwa-mem") {

        // appropriately tagged interleaved FASTQ reads are mapped to the reference
        BWAMEM1_MEM ( BAM2FASTQ.out.reads, bwa, false )
        ch_versions = ch_versions.mix(BWAMEM1_MEM.out.versions)
        aligned_bam = BWAMEM1_MEM.out.bam
    } else {

        // appropriately tagged interleaved FASTQ reads are mapped to the reference
        BWAMEM2_MEM ( BAM2FASTQ.out.reads, bwa, false )
        ch_versions = ch_versions.mix(BWAMEM2_MEM.out.versions)
        aligned_bam = BWAMEM2_MEM.out.bam
    }

    // samblaster is used in order to tag mates information in the BAM file
    // this is used in order to group reads by UMI
    SAMBLASTER ( aligned_bam )
    ch_versions = ch_versions.mix(SAMBLASTER.out.versions)

    // appropriately tagged reads are now grouped by UMI information
    GROUPREADSBYUMI ( SAMBLASTER.out.bam, groupreadsbyumi_strategy )
    ch_versions = ch_versions.mix(GROUPREADSBYUMI.out.versions)

    // using the above created groups, a consensus across reads in the same grou
    // can be called
    // this will emit a consensus BAM file
    CALLUMICONSENSUS ( GROUPREADSBYUMI.out.bam )
    ch_versions = ch_versions.mix(CALLUMICONSENSUS.out.versions)

    emit:
    umibam         = FASTQTOBAM.out.umibam          // channel: [ val(meta), [ bam ] ]
    groupbam       = GROUPREADSBYUMI.out.bam        // channel: [ val(meta), [ bam ] ]
    consensusbam   = CALLUMICONSENSUS.out.bam       // channel: [ val(meta), [ bam ] ]
    versions       = ch_versions                    // channel: [ versions.yml ]
}
