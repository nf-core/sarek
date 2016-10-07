[![Stories in Ready](https://badge.waffle.io/SciLifeLab/CAW.png?label=ready&title=Ready)](https://waffle.io/SciLifeLab/CAW)
# CAW
Nextflow Cancer Analysis Workflow Prototype developed at @SciLifeLab

## Version
0.0.35

## Authors
- Sebastian DiLorenzo (@Sebastian-D)
- Jesper Eisfeldt (@J35P312)
- Maxime Garcia (@MaxUlysse)
- Szilveszter Juhos (@szilvajuhos)
- Max Käller (@gulfshores)
- Malin Larsson (@malinlarsson)
- Björn Nystedt (@bjornnystedt)
- Pall Olason (@pallolason)
- Pelin Sahlén (@pelinakan)

## Quick start
Install [nextflow](http://www.nextflow.io/) and then download the project
```bash
	git clone https://github.com/SciLifeLab/CAW.git
```
Then it's possible to run the full workflow on a data set specified in the tsv file. For example:
```bash
	nextflow run MultiFQtoVC.nf -c milou.config --sample sample.tsv
```
will run the workflow on a small testing dataset. To try on your own data, you just have to edit your own tsv file.

## Usage
```bash
	nextflow run MultiFQtoVC.nf -c <file.config> --sample <file.tsv> --intervals <file.list> [--steps STEP[,STEP]]
```
All variables and parameters are specified in the config and the sample files.

### Steps
To configure which processes will be runned or skipped in the workflow. Different steps to be separated by commas.
Possible values are:
- preprocessing (default, will start workflow with FASTQ files)
- recalibrate (will start workflow with non-recalibrated BAM files)
- skipPreprocessing (will start workflow with recalibrated BAM files)
- MuTect1 (use MuTect1 for VC)
- MuTect2 (use MuTect2 for VC)
- VarDict (use VarDict for VC)
- Strelka (use Strelka for VC)
- HaplotypeCaller (use HaplotypeCaller for normal bams VC)
- Manta (use Manta for SV)
- ascat (use ascat for CNV)

## Nextflow processes
Several processes are run within Nextflow
We divide them for the moment into 2 main steps

### Preprocessing:
- Mapping - Map reads with BWA
- MergeBam - Merge BAMs if multilane samples
- RenameSingleBam - Rename BAM if non-multilane sample
- MarkDuplicates - using Picard
- CreateIntervals - using GATK
- Realign - using GATK
- CreateRecalibrationTable - using GATK
- RecalibrateBam - using GATK

### Variant Calling:
- RunMutect2 - using MuTect2 shipped in GATK v3.6
- concatFiles - merge MuTect2 results
- VarDict - run VarDict on multiple intervals
- VarDictCollatedVCF - merge Vardict results
- RunStrelka - using Strelka 1.0.15
- Manta - run Manta 1.0.0

## TSV file for sample
It's a Tab Separated Value file, based on: "subject status sample lane fastq1 fastq2", "subject status sample bam bai recal" or "subject status sample bam bai"
Quite straight-forward: 
- Subject is the ID of the Patient
- Status is 0 for Normal and 1 for Tumor
- Sample is the Sample ID (It is possible to have more than one tumor sample for each patient)
- Lane is used when the sample is multiplexed on several lanes
- fastq1 is the path to the first pair of the fastq file
- fastq2 is the path to the second pair of the fastq file
- bam is the realigned bam file
- bai is the index
- recal is the recalibration table

### Example TSV file for a normal/tumor pair

In this sample for the normal case there are 3 read groups, and 2 for the tumor. It is recommended to add the absolute path of the paired 
FASTQ files, but relative path should work also. Note, the delimiter is the tab (\t) character:

	G15511	0	normal	normal_1	/samples/G15511/BL/C09DFACXX111207.1_1.fastq.gz	/samples/G15511/BL/C09DFACXX111207.1_2.fastq.gz
	G15511	0	normal	normal_2	/samples/G15511/BL/C09DFACXX111207.2_1.fastq.gz	/samples/G15511/BL/C09DFACXX111207.2_2.fastq.gz
	G15511	0	normal	normal_3	/samples/G15511/BL/C09DFACXX111207.3_1.fastq.gz	/samples/G15511/BL/C09DFACXX111207.3_2.fastq.gz
	G15511	1	tumor	tumor_1	/samples/G15511/T/D0ENMACXX111207.1_1.fastq.gz	/samples/G15511/T/D0ENMACXX111207.1_2.fastq.gz
	G15511	1	tumor	tumor_2	/samples/G15511/T/D0ENMACXX111207.2_1.fastq.gz	/samples/G15511/T/D0ENMACXX111207.2_2.fastq.gz

On the other hand, if you have pre-processed but not-recalibrated BAMs (that is the de-duplicated and realigned BAMs), their indexes and recalibration tables, you should use a structure like:

	G15511	0	normal	/preprocessed/G15511.normal__1.md.real.bam	/preprocessed/G15511.normal__1.md.real.bai	/preprocessed/G15511.normal__1.recal.table
	G15511	1	tumor	/preprocessed/G15511.tumor__1.md.real.bam	/preprocessed/G15511.tumor__1.md.real.bai	/preprocessed/G15511.tumor__1.recal.table

All the files are in the Preprocessing/CreateRecalibrationTable/ directory, and by default a corresponding TSV file is also deposited there. Generally, 
to get MuTect1 and Strelka calls on the preprocessed files should be done by:
```bash
	nextflow -c local.config run MultiFQtoVC.nf --sample Preprocessing/CreateRecalibrationTable_mysample.tsv --steps recalibrate,MuTect1,Strelka
```

The same way, if you have recalibrated BAMs (that is the de-duplicated and realigned and recalibrated BAMs) and their indexes, you should use a structure like:

	G15511	0	normal	/preprocessed/G15511.normal__1.md.real.bam /preprocessed/G15511.normal__1.md.real.bai
	G15511  1	tumor	/preprocessed/G15511.tumor__1.md.real.bam /preprocessed/G15511.tumor__1.md.real.bai

All the files are in he Preprocessing/RecalibrateBam/ directory, and by default a corresponding TSV file is also deposited there. Generally, 
to get MuTect1 and Strelka calls on the recalibrated files should be done by:

```bash
	nextflow -c local.config run MultiFQtoVC.nf --sample Preprocessing/RecalibrateBam_mysample.tsv --steps skipPreprocessing,MuTect1,Strelka
```

## Tools and dependencies
- nextflow 0.17.3
- bwa 0.7.8
- samtools 1.3
- picard 1.118
- GATK 3.6
- R 3.2.3
- gcc 4.9.2
- perl 5.18.4
- strelka 1.0.15
- manta 1.0.0

## Configuration file
We're developing this workflow on UPPMAX cluster milou. So our config file is based on the config file from [NGI RNAseq Nextflow pipeline](https://github.com/SciLifeLab/NGI-RNAseq)

## Use cases

The workflow has three pre-processing options: ```preprocessing```, ```recalibrate``` and skipPreprocessing```. Using
the ```preprocessing``` directive one will have a pair of mapped, deduplicated and recalibrated BAM files in the
```Preprocessing/RecalibrateBam/``` directory. Furthermore, during this process a deduplicated and realigned BAM file is
created alongside with its recalibration table in the ```Preprocessing/CreateRecalibrationTable/``` directory. This
latter BAM is usually much smaller and easier to backup together with its recalibration table. Making recalibrated BAMs
from this input data is only one step. This is the usual option you have to give when you are starting from raw FASTQ
data:

```
nextflow run MultiFQtoVC.nf -c milou.config --sample mysample.tsv
```

Preprocessing will start by default, you do not have to give any additional steps, only the TSV file describing the
sample (see below). 

In the provided milou.config file we are defining the intervals file as well, this is used to
define regions for variant call and realignment (in a scatter and gather fashion when possible). The intervals are
chromosomes cut at their centromeres (so each chromosome arm processed separately) also additional unassigned contigs.
We are ignoring the hs37d5 contig that contains concatenated decoy sequences. 

### Starting from raw FASTQ - having pair of FASTQ files for normal and tumor samples (one lane for each sample)

Ide teszunk valami ertelmes szoveget majd

### Starting from raw FASTQ - having multiple lanes (reads groups)

Ide meg a szoveg tobbi resze jon majd

### Starting from a preprocessed BAM file with recalibration table

Ide majd hogy itt mit kell csinalni

### Starting from a recalibrated BAM file

At this step we are assuming that all the required preprocessing steps (alignment, deduplication, ..., recalibration) is over, we only want to run variant callers or other tools using recalibrated BAMs.


## Test data
For quick, easy testing of the workflow, we need a small, but representative
data set as input. We went for a full depth coverage, small chunk of the genome. "Seba" suggested using [TCGA benchmarking data][TCGA], which is already downloaded to Uppmax.

Grabbing 100kb from chr1, both tumor and normal, sorting by read name:
```bash
	samtools view -bh /sw/data/uppnex/ToolBox/TCGAbenchmark/f0eaa94b-f622-49b9-8eac-e4eac6762598/G15511.HCC1143_BL.1.bam 1:100000-200000 | samtools sort -n -m 4G - tcga.cl.normal.part.nsort.bam
	samtools view -bh /sw/data/uppnex/ToolBox/TCGAbenchmark/ad3d4757-f358-40a3-9d92-742463a95e88/G15511.HCC1143.1.bam 1:100000-200000 | samtools sort -n -m 4G - tcga.cl.tumor.part.nsort.bam
```

Then convert to fastq:
```bash
	bedtools bamtofastq -i tcga.cl.normal.part.nsort.bam -fq tcga.cl.normal.part.nsort_R1.fastq -fq2 tcga.cl.normal.part.nsort_R2.fastq
	bedtools bamtofastq -i tcga.cl.tumor.part.nsort.bam -fq tcga.cl.tumor.part.nsort_R1.fastq -fq2 tcga.cl.tumor.part.nsort_R2.fastq
```
And finally gzipping...

[TCGA]: https://cghub.ucsc.edu/datasets/benchmark_download.html
