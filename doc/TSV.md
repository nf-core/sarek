# TSV file for sample
It's a Tab Separated Value file, based on: `subject gender status sample lane fastq1 fastq2` or `subject gender status sample bam bai`
Quite straight-forward:
- `subject` is the ID of the Patient
- `gender` is the gender of the Patient, (XX or XY)
- `status` is the status of the Patient, (0 for Normal or 1 for Tumor)
- `sample` is the Sample ID (It is possible to have more than one tumor sample for each patient)
- `lane` is used when the sample is multiplexed on several lanes
- `fastq1` is the path to the first pair of the fastq file
- `fastq2` is the path to the second pair of the fastq file
- `bam` is the bam file
- `bai` is the index

# Example TSV file for a normal/tumor pair with FASTQ files

In this sample for the normal case there are 3 read groups, and 2 for the tumor. It is recommended to add the absolute path of the paired 
FASTQ files, but relative path should work also. Note, the delimiter is the tab (\t) character:
```
G15511	XX	0	normal	normal_1	pathToFiles/C09DFACXX111207.1_1.fastq.gz	pathToFiles/C09DFACXX111207.1_2.fastq.gz
G15511	XX	0	normal	normal_2	pathToFiles/C09DFACXX111207.2_1.fastq.gz	pathToFiles/C09DFACXX111207.2_2.fastq.gz
G15511	XX	0	normal	normal_3	pathToFiles/C09DFACXX111207.3_1.fastq.gz	pathToFiles/C09DFACXX111207.3_2.fastq.gz
G15511	XX	1	tumor	tumor_1	pathToFiles/D0ENMACXX111207.1_1.fastq.gz	pathToFiles/D0ENMACXX111207.1_2.fastq.gz
G15511	XX	1	tumor	tumor_2	pathToFiles/D0ENMACXX111207.2_1.fastq.gz	pathToFiles/D0ENMACXX111207.2_2.fastq.gz
```

# Example TSV file for a normal/tumor pair with BAM files
On the other hand, if you have BAMs (T/N pairs that were not realigned together) and their indexes, you should use a structure like:
```
G15511	XX	0	normal	pathToFiles/G15511.normal__1.md.real.bam	pathToFiles/G15511.normal__1.md.real.bai
G15511	XX	1	tumor	pathToFiles/G15511.tumor__1.md.real.bam pathToFiles/G15511.tumor__1.md.real.bai
```
All the files will be created in the Preprocessing/NonRealigned/ directory, and by default a corresponding TSV file will also be deposited there. Generally, getting MuTect1 and Strelka calls on the preprocessed files should be done by:
```bash
nextflow run SciLifeLab/CAW --sample Preprocessing/NonRealigned/mysample.tsv --steps realign,MuTect1,Strelka
```
# Example TSV file for a normal/tumor pair with recalibrated BAM files
The same way, if you have recalibrated BAMs (T/N pairs that were realigned together) and their indexes, you should use a structure like:
```
G15511	XX	0	normal	pathToFiles/G15511.normal__1.md.real.bam	pathToFiles/G15511.normal__1.md.real.bai
G15511	XX	1	tumor	pathToFiles/G15511.tumor__1.md.real.bam	pathToFiles/G15511.tumor__1.md.real.bai
```
All the files will be in he Preprocessing/Recalibrated/ directory, and by default a corresponding TSV file will also be deposited there. Generally, getting MuTect1 and Strelka calls on the recalibrated files should be done by:

```bash
nextflow run SciLifeLab/CAW --sample Preprocessing/Recalibrated/mysample.tsv --steps skipPreprocessing,MuTect1,Strelka
```
