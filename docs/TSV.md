# TSV file for sample

It's a Tab Separated Value file, based on: `subject gender status sample lane fastq1 fastq2` or `subject gender status sample bam bai` Quite straight-forward:

- `subject` designate the subject, it should be the ID of the Patient, or if you don't have one, il could be the Normal ID Sample.
- `gender` is the gender of the Patient, (XX or XY)
- `status` is the status of the Patient, (0 for Normal or 1 for Tumor)
- `sample` designate the Sample, it should be the ID of the Sample (it is possible to have more than one tumor sample for each patient)
- `lane` is used when the sample is multiplexed on several lanes
- `fastq1` is the path to the first pair of the fastq file
- `fastq2` is the path to the second pair of the fastq file
- `bam` is the bam file
- `bai` is the index

All examples are given for a normal/tumor pair. If no tumors are listed in the TSV file, then the workflow will proceed as if it was a single normal sample instead of a normal/tumor pair.

# Example TSV file for a normal/tumor pair with FASTQ files

In this sample for the normal case there are 3 read groups, and 2 for the tumor. It is recommended to add the absolute path of the paired FASTQ files, but relative path should work also. Note, the delimiter is the tab (\t) character:

```
G15511    XX    0    C09DFN    C09DF_1    pathToFiles/C09DFACXX111207.1_1.fastq.gz    pathToFiles/C09DFACXX111207.1_2.fastq.gz
G15511    XX    0    C09DFN    C09DF_2    pathToFiles/C09DFACXX111207.2_1.fastq.gz    pathToFiles/C09DFACXX111207.2_2.fastq.gz
G15511    XX    0    C09DFN    C09DF_3    pathToFiles/C09DFACXX111207.3_1.fastq.gz    pathToFiles/C09DFACXX111207.3_2.fastq.gz
G15511    XX    1    D0ENMT    D0ENM_1    pathToFiles/D0ENMACXX111207.1_1.fastq.gz    pathToFiles/D0ENMACXX111207.1_2.fastq.gz
G15511    XX    1    D0ENMT    D0ENM_2    pathToFiles/D0ENMACXX111207.2_1.fastq.gz    pathToFiles/D0ENMACXX111207.2_2.fastq.gz
```

# Example TSV file for a normal/tumor pair with BAM files

On the other hand, if you have BAMs (T/N pairs that were not realigned together) and their indexes, you should use a structure like:

```
G15511    XX    0    C09DFN    pathToFiles/G15511.C09DFN.md.real.bam    pathToFiles/G15511.C09DFN.md.real.bai
G15511    XX    1    D0ENMT    pathToFiles/G15511.D0ENMT.md.real.bam pathToFiles/G15511.D0ENMT.md.real.bai
```

All the files will be created in the Preprocessing/NonRealigned/ directory, and by default a corresponding TSV file will also be deposited there. Generally, getting MuTect1 and Strelka calls on the preprocessed files should be done by:

```bash
nextflow run SciLifeLab/Sarek/main.nf --sample Preprocessing/NonRealigned/mysample.tsv --step realign
nextflow run SciLifeLab/Sarek/somaticVC.nf --sample Preprocessing/Recalibrated/mysample.tsv --tools Mutect2,Strelka
```

# Example TSV file for a normal/tumor pair with recalibrated BAM files

The same way, if you have recalibrated BAMs (T/N pairs that were realigned together) and their indexes, you should use a structure like:

```
G15511    XX    0    C09DFN    pathToFiles/G15511.C09DFN.md.real.bam    pathToFiles/G15511.C09DFN.md.real.bai
G15511    XX    1    D0ENMT    pathToFiles/G15511.D0ENMT.md.real.bam    pathToFiles/G15511.D0ENMT.md.real.bai
```

All the files will be in he Preprocessing/Recalibrated/ directory, and by default a corresponding TSV file will also be deposited there. Generally, getting MuTect1 and Strelka calls on the recalibrated files should be done by:

```bash
nextflow run SciLifeLab/Sarek/somaticVC.nf --sample Preprocessing/Recalibrated/mysample.tsv --tools Mutect2,Strelka
```

--------------------------------------------------------------------------------

[![](images/SciLifeLab_logo.png "SciLifeLab")][scilifelab-link]
[![](images/NGI_logo.png "NGI")][ngi-link]
[![](images/NBIS_logo.png "NBIS")][nbis-link]

[nbis-link]: https://www.nbis.se/
[ngi-link]: https://ngisweden.scilifelab.se/
[scilifelab-link]: https://www.scilifelab.se/
