# NF-CORE fastq_align_dna

```mermaid
flowchart LR
reads & index & aligner --> sort{Sort?} --> aligners
aligners --> bam & bai & reports & versions

subgraph aligners
  bowtie2
  bwamem
  bwamem2
  dragmap
  snap
end

```
