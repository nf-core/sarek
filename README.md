# CAW
Nextflow Cancer Analysis Workflow Prototype developed at @SciLifeLab

## Version
0.0.3

## Authors
- Pelin Akan (@pelinakan)
- Jesper Eisfeldt (@J35P312)
- Maxime Garcia (@MaxUlysse)
- Szilveszter Juhos (@szilvajuhos)
- Max Käller
- Malin Larsson (@malinlarsson)
- Björn Nystedt (@bjornnystedt)
- Pall Olason (@pallolason)

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
- nopreprocessing (will start workflow with recalibrated BAM files)
- MuTect2 (use MuTect2 for VC)
- VarDict (use VarDict for VC)
- Strelka (use Strelka for VC)
- Manta (use Manta for SV)

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
- Manta - run Manta 0.27.1

## TSV file for sample
It's a Tab Separated Value file, based on: "subject status sample lane fastq1 fastq2" or "subject status sample bam bai recal"
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

## tools and dependencies
nextflow 0.17.3
bwa 0.7.8
samtools 1.3
picard 1.118
java sun_jdk 1.8.0_92
GATK 3.6
R 3.2.3
gcc 4.9.2
java sun_jdk1.8.0_40
perl 5.18.4
strelka 1.0.15
manta 0.27.1

## Configuration file
We're developing this workflow on UPPMAX cluster milou. So our config file is based on the config file from [NGI RNAseq Nextflow pipeline](https://github.com/SciLifeLab/NGI-RNAseq)

## Test data
For quick, easy testing of the workflow, we need a small, but representative
data set as input. We went for a full depth coverage, small chunk of the genome. "Seba" suggested using [TCGA benchmarking data][TCGA], which is already downloaded to Uppmax.

Grabbing 100kb from chr1, both tumor and normal, sorting by read name:
```
$ samtools view -bh /sw/data/uppnex/ToolBox/TCGAbenchmark/f0eaa94b-f622-49b9-8eac-e4eac6762598/G15511.HCC1143_BL.1.bam 1:100000-200000 | samtools sort -n -m 4G - tcga.cl.normal.part.nsort.bam
$ samtools view -bh /sw/data/uppnex/ToolBox/TCGAbenchmark/ad3d4757-f358-40a3-9d92-742463a95e88/G15511.HCC1143.1.bam 1:100000-200000 | samtools sort -n -m 4G - tcga.cl.tumor.part.nsort.bam
```

Then convert to fastq:
```
$ bedtools bamtofastq -i tcga.cl.normal.part.nsort.bam -fq tcga.cl.normal.part.nsort_R1.fastq -fq2 tcga.cl.normal.part.nsort_R2.fastq
$ bedtools bamtofastq -i tcga.cl.tumor.part.nsort.bam -fq tcga.cl.tumor.part.nsort_R1.fastq -fq2 tcga.cl.tumor.part.nsort_R2.fastq
```
And finally gzipping...

[TCGA]: https://cghub.ucsc.edu/datasets/benchmark_download.html
