# TSV file for sample(s)

Input files for Sarek can be specified using a TSV file given to the `--sample` command.
The TSV file is a Tab Separated Value file with columns:
* `subject gender status sample lane fastq1 fastq2` for step `mapping` with paired-end FASTQs
* `subject gender status sample lane bam` for step `mapping` with unmapped BAMs
* `subject gender status sample bam bai recaltable` for step `recalibrate` with BAMs
* `subject gender status sample bam bai` for step `variantcalling` with BAMs

The content of these columns is quite straight-forward:
* `subject` designate the subject, it should be the ID of the Patient
* `gender` is the gender of the Patient, (XX or XY)
* `status` is the status of the Patient, (0 for Normal or 1 for Tumor)
* `sample` designate the Sample, it should be the ID of the Sample (it is possible to have more than one tumor sample for each patient)
* `lane` is used when the sample is multiplexed on several lanes
* `fastq1` is the path to the first pair of the fastq file
* `fastq2` is the path to the second pair of the fastq file
* `bam` is the bam file
* `bai` is the bam index file
* `recaltable` is the recalibration table

It is recommended to add the absolute path of the files, but relative path should work also.
Note, the delimiter is the tab (`\t`) character:

All examples are given for a normal/tumor pair.
If no tumors are listed in the TSV file, then the workflow will proceed as if it is a normal sample instead of a normal/tumor pair.

Sarek will output results is a different directory for each sample.
If multiple samples are specified in the TSV file, Sarek will consider all files to be from different samples.
Multiple TSV files can be specified if the path is enclosed in quotes.

Somatic variant calling output will be in a specific directory for each normal/tumor pair.

# Example TSV file for a normal/tumor pair with FASTQ files

In this sample for the normal case there are 3 read groups, and 2 for the tumor.

```
G15511    XX    0    C09DFN    C09DF_1    pathToFiles/C09DFACXX111207.1_1.fastq.gz    pathToFiles/C09DFACXX111207.1_2.fastq.gz
G15511    XX    0    C09DFN    C09DF_2    pathToFiles/C09DFACXX111207.2_1.fastq.gz    pathToFiles/C09DFACXX111207.2_2.fastq.gz
G15511    XX    0    C09DFN    C09DF_3    pathToFiles/C09DFACXX111207.3_1.fastq.gz    pathToFiles/C09DFACXX111207.3_2.fastq.gz
G15511    XX    1    D0ENMT    D0ENM_1    pathToFiles/D0ENMACXX111207.1_1.fastq.gz    pathToFiles/D0ENMACXX111207.1_2.fastq.gz
G15511    XX    1    D0ENMT    D0ENM_2    pathToFiles/D0ENMACXX111207.2_1.fastq.gz    pathToFiles/D0ENMACXX111207.2_2.fastq.gz
```

# Example TSV file for a normal/tumor pair with BAM files

In this sample for the normal case there are 3 read groups, and 2 for the tumor.

```
G15511    XX    0    C09DFN    C09DF_1    pathToFiles/C09DFAC_1.bam
G15511    XX    0    C09DFN    C09DF_2    pathToFiles/C09DFAC_2.bam
G15511    XX    0    C09DFN    C09DF_3    pathToFiles/C09DFAC_3.bam
G15511    XX    1    D0ENMT    D0ENM_1    pathToFiles/D0ENMAC_1.bam
G15511    XX    1    D0ENMT    D0ENM_2    pathToFiles/D0ENMAC_2.bam
```

# Example TSV file for a normal/tumor pair with recalibrated BAM files

The same way, if you have recalibrated BAMs and their indexes, you should use a structure like:

```
G15511    XX    0    C09DFN    pathToFiles/G15511.C09DFN.md.real.bam    pathToFiles/G15511.C09DFN.md.real.bai
G15511    XX    1    D0ENMT    pathToFiles/G15511.D0ENMT.md.real.bam    pathToFiles/G15511.D0ENMT.md.real.bai
```

## Input FASTQ file name best practices

The input folder, containing the FASTQ files for one individual (ID) should be organized into one subfolder for every sample.
All fastq files for that sample should be collected here.

```
ID
+--sample1
+------sample1_lib_flowcell-index_lane_R1_1000.fastq.gz
+------sample1_lib_flowcell-index_lane_R2_1000.fastq.gz
+------sample1_lib_flowcell-index_lane_R1_1000.fastq.gz
+------sample1_lib_flowcell-index_lane_R2_1000.fastq.gz
+--sample2
+------sample2_lib_flowcell-index_lane_R1_1000.fastq.gz
+------sample2_lib_flowcell-index_lane_R2_1000.fastq.gz
+--sample3
+------sample3_lib_flowcell-index_lane_R1_1000.fastq.gz
+------sample3_lib_flowcell-index_lane_R2_1000.fastq.gz
+------sample3_lib_flowcell-index_lane_R1_1000.fastq.gz
+------sample3_lib_flowcell-index_lane_R2_1000.fastq.gz
```

Fastq filename structure:

- `sample_lib_flowcell-index_lane_R1_1000.fastq.gz` and
- `sample_lib_flowcell-index_lane_R2_1000.fastq.gz`

Where:

- `sample` = sample id
- `lib` = indentifier of libaray preparation
- `flowcell` = identifyer of flow cell for the sequencing run
- `lane` = identifier of the lane of the sequencing run

Read group information will be parsed from fastq file names according to this:

- `RGID` = "sample_lib_flowcell_index_lane"
- `RGPL` = "Illumina"
- `PU` = sample
- `RGLB` = lib

# Path to a FASTQ directory for a single normal sample

Input files for Sarek can be specified using the path to a FASTQ directory given to the `--sample` command only with the `mapping` step.

# VCF files for annotation

Input files for Sarek can be specified using the path to a VCF directory given to the `--sample` command only with the `annotate` step.
Multiple VCF files can be specified if the path is enclosed in quotes.
