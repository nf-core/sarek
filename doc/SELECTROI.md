# Selecting region of interest (ROI) from BAM files

## Purpose
It is difficult to look at the alignments of huge WGS files, and most of the data is usually have little interest. The purpose of the 'selectROI.py' 
script is to cut out only the interesting parts from the large file, and save it into a smaller file (or files) that can be stored on and visualize
even on a laptop or desktop machine. 

## Prerequisites

The script was tested with python 2.7 but should work with python 3.x . You have to have [samtools](http://www.htslib.org/) installed or added by a module system. 
The BAM alignment files have to be indexed and the BAM index file have to be next to the alignment file in the same directory.
My experience was that selecting 100K variants and 1200 genes from a 80G BAM file resulted in a 8G file, that can be managed on a laptop.
Some parts of the code were stolen/reused from the [interlap](https://github.com/brentp/interlap) python module.

## Command-line parameters

 * **-a alignment.bam** Input alignment file(s): a single indexed BAM, or a list of indexed BAM files separated by commas. This is a mandatory argument,
 at least one input alignment BAM file is needed.

 * **-b BEDfile.bed** A single [BED file](https://genome.ucsc.edu/FAQ/FAQformat.html#format1) or a set of BED files separated by commas. Usually this is a list
 of interesting genes or regions. Reads between the starting and the ending position of the listed genes are exported. Duplicates and overlaps are handled gracefully,
 meaning that duplicates and overlaps are collapsed into a single region. This is an optional argument, but at least a VCF or a BED file have to be specified 
 in the command line by either `-v` or `-b`.

 * **-h help** Short command-line help

 * **-m maxRegions** Max number of regions processed in a single chunk. Optional, default is 150. This script is using `samtools`, and forks a new
 process if there are at least this number of regions to be passed into the command line. If you have less memory, decrease this value. If you are 
 increasing this value too much, the underlying shell will bail out and complaining about a too long command-line argument list.

 * **-p prefix** Prefix of the output file(s). Optional, default is **ROI.**. Input file `testfile.bam` will produce an output `ROI.testfile.bam` .
 
 * **-r readlength** Padding, the width (see the `-w` option) is expanded with this pad that is corresponding to read length. Optional, default value
 is 150. With a 200 width and 150 padding an 700 bps wide region is exported around a single SNP coordinate. Regions that are closer to each other than
 this pad value are merged. I.e. if there are two SNPs 500 bps distance next to each other, using the default values there would be a 100 bps gap in
 the output. However, because their regions are expanded with the pad parameter, they will be merged, and all reads between the two (and additional
 reads preceding and following their coordinates) are exported. Using default values two SNPs have to be 200+150+150+200 +1  = 701 bps apart to be 
 exported into two separate regions. This also means that there will be a coverage dip between the two variants. 

 * **-t threads** Number of CPU threads/cores to use, optional (default is 8)

 * **-v VCFfile.vcf** A single [vcf file](https://samtools.github.io/hts-specs/) or a set of VCF files separated by commas. It is expected that variants in 
 this file are are important for some reason, as for each VCF line the chromosome and start coordinate is extracted, and all the reads in the surrounding of
 the variant are exported to the output BAM file. The breadth of the selection is controlled by the **-w** and **-p** parameters. When supplying more than one
 VCFs, use commas like `-v file1.vcf,file2.vcf,/somewhere/else/file3.vcf` . Duplicate coordinates are filtered out, so reads around close/overlapping variants
 are exported only once. This is an optional argument, but at least a VCF or a BED file have to be specified in the command line by either `-v` or `-b`.

 * **-w width** Region width/breadth for VCF records. Optional, default is 200. For each VCF record reads starting/ending this number of bases before 
 and after the coordinate of the record respectively are exported. Actually this values is expanded with the **pad** parameter (see the `-r` option). 


## Examples
Be sure you can access `samtools` either by installing it or load the module.

 - Extracting all the genes listed in a BED from a BAM file:

```bash
python selectROI.py -a alignment.bam -b genes.bed
```

 - Ditto, but extracting all the genes listed in a BED from a list of BAM files:

```bash
python selectROI.py -a alignment1.bam,alignment2.bam,alignment3.bam -b genes.bed
```

 - Extracting regions around germline and somatic variants and genes from a BAM file:

```bash
python selectROI.py -a alignment.bam -b genes.bed -v germline.vcf,somatic.vcf 
```

 - Getting reads for a normal/tumour pair from a list of VCFs and genes:

```bash
python selectROI.py -a normal.bam,tumour.bam -b genes.bed -v callset1.vcf,callset2.vcf,callset3.vcfq
```
