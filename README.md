[![Stories in Ready](https://badge.waffle.io/SciLifeLab/CAW.png?label=ready&title=Ready)](https://waffle.io/SciLifeLab/CAW)

# CAW
Nextflow Cancer Analysis Workflow Prototype developed at the [National Genomics Infastructure](https://ngisweden.scilifelab.se/)
at [SciLifeLab Stockholm](https://www.scilifelab.se/platforms/ngi/), Sweden.

## Version
0.8.2

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

##Installation

###Install Nextflow
To use this pipeline, you need to have a working version of NextFlow installed. You can find more information about this pipeline tool at [nextflow.io](http://www.nextflow.io/). The typical installation
of NextFlow looks like this:
```bash
curl -fsSL get.nextflow.io | bash
mv ./nextflow ~/bin
```

#### UPPMAX
If you're running on a Swedish UPPMAX cluster you can load NextFlow as an environment module instead:
```bash
module load Nextflow
```
The first time you load this you will get a warning about setting environment variables. To automatically set these at login, you can add the following lines to your `~/.bashrc` file:
```bash
export NXF_LAUNCHBASE=$SNIC_TMP
export NXF_TEMP=$SNIC_TMP
export NXF_OPTS='-Xms1g -Xmx4g '
```

### NextFlow configuration
Next, you need to set up a config file so that NextFlow knows how to run and where to find reference
indexes. You can find an example configuration file for UPPMAX (milou) with this repository:
[`example_uppmax_config`](https://github.com/SciLifeLab/NGI-RNAseq/blob/master/example_uppmax_config).

Copy this file to `~/.nextflow/config` and edit the line `'-A b2013064'` to contain your own UPPMAX project
identifier instead.

It is entirely possible to run this pipeline on other clusters - just note that you may need to customise
the `process` environment (eg. if you're using a cluster system other than SLURM) and the paths to reference
files.

### Workflow installation
This workflow itself needs no installation - NextFlow will automatically fetch it from GitHub when run if `SciLifeLab/CAW` is specified as the workflow name.

If you prefer, you can download the files yourself from GitHub and run them directly:
```bash
git clone https://github.com/SciLifeLab/CAW
nextflow run CAW/main.nf
```

It's possible to run the full workflow on a data set specified in the tsv file. For example:
```bash
nextflow run SciLifeLab/CAW -c milou.config --sample sample.tsv
```
will run the workflow on a small testing dataset. To try on your own data, you just have to edit your own tsv file.

## Usage
```bash
nextflow run SciLifeLab/CAW -c <file.config> --sample <file.tsv> --intervals <file.list> [--steps STEP[,STEP]]
```
All variables and parameters are specified in the config and the sample files.

## Cleaning
Use nextflow clean -f to remove everything contained in the work directories. Do not worry, non-recalibrated and recalibrated bams and index are stored in the `Preprocessing` directory. And variant calling files are stored in the `VariantCalling` directory.
```bash
nextflow clean -f
```

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
It's a Tab Separated Value file, based on: `subject status sample lane fastq1 fastq2`, `subject status sample bam bai recal` or `subject status sample bam bai`
Quite straight-forward:
- `subject` is the ID of the Patient
- `status` is the status of the Patient, (0 for Normal and 1 for Tumor)
- `sample` is the Sample ID (It is possible to have more than one tumor sample for each patient)
- `lane` is used when the sample is multiplexed on several lanes
- `fastq1` is the path to the first pair of the fastq file
- `fastq2` is the path to the second pair of the fastq file
- `bam` is the realigned bam file
- `bai` is the index
- `recal` is the recalibration table

### Example TSV file for a normal/tumor pair

In this sample for the normal case there are 3 read groups, and 2 for the tumor. It is recommended to add the absolute path of the paired 
FASTQ files, but relative path should work also. Note, the delimiter is the tab (\t) character:
```
G15511	0	normal	normal_1	/samples/G15511/BL/C09DFACXX111207.1_1.fastq.gz	/samples/G15511/BL/C09DFACXX111207.1_2.fastq.gz
G15511	0	normal	normal_2	/samples/G15511/BL/C09DFACXX111207.2_1.fastq.gz	/samples/G15511/BL/C09DFACXX111207.2_2.fastq.gz
G15511	0	normal	normal_3	/samples/G15511/BL/C09DFACXX111207.3_1.fastq.gz	/samples/G15511/BL/C09DFACXX111207.3_2.fastq.gz
G15511	1	tumor	tumor_1	/samples/G15511/T/D0ENMACXX111207.1_1.fastq.gz	/samples/G15511/T/D0ENMACXX111207.1_2.fastq.gz
G15511	1	tumor	tumor_2	/samples/G15511/T/D0ENMACXX111207.2_1.fastq.gz	/samples/G15511/T/D0ENMACXX111207.2_2.fastq.gz
```
On the other hand, if you have pre-processed but not-recalibrated BAMs (that is the de-duplicated and realigned BAMs), their indexes and recalibration tables, you should use a structure like:
```
G15511	0	normal	pathToFiles/G15511.normal__1.md.real.bam	pathToFiles/G15511.normal__1.md.real.bai	pathToFiles/G15511.normal__1.recal.table
G15511	1	tumor	pathToFiles/G15511.tumor__1.md.real.bam	pathToFiles/G15511.tumor__1.md.real.bai	pathToFiles/G15511.tumor__1.recal.table
```
All the files will be created in the Preprocessing/NonRecalibrated/ directory, and by default a corresponding TSV file will also be deposited there. Generally, getting MuTect1 and Strelka calls on the preprocessed files should be done by:
```bash
nextflow -c local.config run MultiFQtoVC.nf --sample Preprocessing/NonRecalibrated/mysample.tsv --steps recalibrate,MuTect1,Strelka
```

The same way, if you have recalibrated BAMs (that is the de-duplicated and realigned and recalibrated BAMs) and their indexes, you should use a structure like:
```
G15511	0	normal	pathToFiles/G15511.normal__1.md.real.bam pathToFiles/G15511.normal__1.md.real.bai
G15511  1	tumor	pathToFiles/G15511.tumor__1.md.real.bam pathToFiles/G15511.tumor__1.md.real.bai
```
All the files will be in he Preprocessing/Recalibrated/ directory, and by default a corresponding TSV file will also be deposited there. Generally, getting MuTect1 and Strelka calls on the recalibrated files should be done by:

```bash
nextflow -c local.config run MultiFQtoVC.nf --sample Preprocessing/Recalibrated/mysample.tsv --steps skipPreprocessing,MuTect1,Strelka
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

The workflow has three pre-processing options: `preprocessing`, `recalibrate` and `skipPreprocessing`. Using the `preprocessing` directive one will have a pair of mapped, deduplicated and recalibrated BAM files in the `Preprocessing/Recalibrated/` directory. Furthermore, during this process a deduplicated and realigned BAM file is created alongside with its recalibration table in the `Preprocessing/NonRecalibrated/` directory. This latter BAM is usually much smaller and easier to backup together with its recalibration table. Making recalibrated BAMs from this input data is only one step. This is the usual option you have to give when you are starting from raw FASTQ data:

```bash
nextflow run MultiFQtoVC.nf -c milou.config --sample mysample.tsv
```

Preprocessing will start by default, you do not have to give any additional steps, only the TSV file describing the sample (see below). 

In the provided milou.config file we are defining the intervals file as well, this is used to define regions for variant call and realignment (in a scatter and gather fashion when possible). The intervals are chromosomes cut at their centromeres (so each chromosome arm processed separately) also additional unassigned contigs. We are ignoring the hs37d5 contig that contains concatenated decoy sequences. 

During processing steps a `trace.txt` and a `timeline.html` file is generated automatically. These files contain statistics about resources used and processes finished. If you start a new flow or restart/resume a sample, the previous version will be renamed as `trace.txt.1` and `timeline.html.1` respectively. Also, older version are renamed with incremented numbers.

### Starting from raw FASTQ - having pair of FASTQ files for normal and tumor samples (one lane for each sample)

The workflow should be started in this case with the smallest set of options as written above:
```bash
nextflow run MultiFQtoVC.nf -c milou.config --sample mysample.tsv
```
The TSV file should have at least two tab-separated lines:
```
SUBJECT_ID	0	normal	1	/samples/normal_1.fastq.gz	/samples/normal_2.fastq.gz
SUBJECT_ID	1	tumor	1	/samples/tumor_1.fastq.gz	/samples/tumor_2.fastq.gz
```
The columns are:
1. Subject id
2. identifier: 0 if normal, 1 if tumor
3. text id: actual text representation of the type of the sample
4. read group ID: it is irrelevant in this simple case, should to be 1
5. first set of reads
6. second set of reads

### Starting from raw FASTQ - having pair of FASTQ files for normal, tumor and relapse samples (one lane for each sample)

The workflow command line is just the same as before, but the TSV contains an extra line. You can see the second column is used to distinguish normal and tumor, there is no extra identifier for relapse. You can add as many relapse samples as many you have, providing their name in the third column is different. Each will be compared to the normal
one-by-one.
```
SUBJECT_ID	0	normal	1	/samples/normal_1.fastq.gz	/samples/normal_2.fastq.gz
SUBJECT_ID	1	tumor	1	/samples/tumor_1.fastq.gz	/samples/tumor_2.fastq.gz
SUBJECT_ID	1	relapse	1	/samples/relapse_1.fastq.gz	/samples/relapse_2.fastq.gz
```
### Starting from raw FASTQ - having multiple lanes (reads groups)

Usually there are more read groups - sequencing lanes - for a single sequencing run, and in a flowcell different lanes have to be recalibrated separately. This is captured in the TSV file only in the following manner, adding read group numbers or IDs in the fourth column. Obviously, if you do not have relapse samples, you can leave out those lines.
```
SUBJECT_ID	0	normal	1	/samples/normal1_1.fastq.gz	/samples/normal1_2.fastq.gz
SUBJECT_ID	0	normal	2	/samples/normal2_1.fastq.gz	/samples/normal2_2.fastq.gz
SUBJECT_ID	1	tumor	3	/samples/tumor3_1.fastq.gz	/samples/tumor3_2.fastq.gz
SUBJECT_ID	1	tumor	4	/samples/tumor4_1.fastq.gz	/samples/tumor4_2.fastq.gz
SUBJECT_ID	1	tumor	5	/samples/tumor5_1.fastq.gz	/samples/tumor5_2.fastq.gz
SUBJECT_ID	1	relapse	7	/samples/relapse7_1.fastq.gz	/samples/relapse7_2.fastq.gz
SUBJECT_ID	1	relapse	9	/samples/relapse9_1.fastq.gz	/samples/relapse9_2.fastq.gz
```

### Starting from recalibration

NGI Production in the previous years delivered many preprocessed samples; these BAM files are not recalibrated, but their recalibration table was also delivered together with the alignments. To have recalibrated BAMs (suitable for
variant calling) is a single additional step:
```bash
nextflow run MultiFQtoVC.nf -c milou.config --sample mysample.tsv --steps recalibrate
```
And the corresponding TSV file should be like:
```
SUBJECT_ID	0	normal	1	/samples/normal.bam	/samples/normal.recal.table
SUBJECT_ID	1	tumor	1	/samples/tumor.bam	/samples/tumor.recal.table
```
At the end of this step you should have recalibrated BAM files in the `Preprocessing/Recalibrated/` directory.

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
