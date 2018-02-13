# Selecting region of interest (ROI) from BAM files

## Purpose
It is difficult to look at the alignments of huge WGS files, and most of the data is usually have little interest. The purpose of the 'selectROI.py' 
script is to cut out only the interesting parts from the large file, and save it into a smaller file (or files) that can be stored on and visualize
even on a laptop or desktop machine. 

## Prerequisites

The script was tested with python 2.7 but should work with python 3.x . You have to have [samtools](http://www.htslib.org/) installed or added by a module system.

My experience was that selecting 100K variants and 1200 genes from a 80G BAM file resulted in a 8G file, that can be managed on a laptop.

Some parts of the code were stolen/reused from the [interlap](https://github.com/brentp/interlap) python module.

## Command-line parameters

## Examples 




```bash
module load bioinfo-tools samtools6.4 bwa
samtools faidx /sw/data/uppnex/ToolBox/ReferenceAssemblies/hg38make/bundle/2.8/b37/human_g1k_v37_decoy.fasta $(< tiny.list) > intervals.fasta
bwa mem /sw/data/uppnex/ToolBox/hg38bundle/Homo_sapiens_assembly38.fasta intervals.fasta > intervals.sam
grep -v '^@' intervals.sam | awk '{printf("%s:%d-%d\n", $3, $4, $4+$6-1)}' > tiny-GRCh38.list
```

--------------------------------------------------------------------------------

[![](images/SciLifeLab_logo.png "SciLifeLab")][scilifelab-link]
[![](images/NGI_logo.png "NGI")][ngi-link]
[![](images/NBIS_logo.png "NBIS")][nbis-link]

[nbis-link]: https://www.nbis.se/
[ngi-link]: https://ngisweden.scilifelab.se/
[scilifelab-link]: https://www.scilifelab.se/
