# Project folder structure

The workflow is started for a sample, or a set of samples from the same person.
Each different physical samples is identified by its own ID.
For example in a Tumor/Normal settings, this ID could correspond to "Normal", "Tumor 1", "Tumor 2" etc corresponding to all physical samples from the same cancer patient.

Below is an overview of the intended folder structure for an analyzed project.

![Project folder structure](images/folder_structure.jpg "Folder structure")

# Input Fastq file name conventions

The input folder, containing the fastq files for one ID (Individual) should be organized into one subfolder for every sample. All fastq files for that sample should be collected here.

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

--------------------------------------------------------------------------------

[![](images/SciLifeLab_logo.png "SciLifeLab")][scilifelab-link]
[![](images/NGI_logo.png "NGI")][ngi-link]
[![](images/NBIS_logo.png "NBIS")][nbis-link]

[nbis-link]: https://www.nbis.se/
[ngi-link]: https://ngisweden.scilifelab.se/
[scilifelab-link]: https://www.scilifelab.se/
