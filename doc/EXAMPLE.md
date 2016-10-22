# Example TSV file for a normal/tumor pair with FASTQ files

In this sample for the normal case there are 3 read groups, and 2 for the tumor. It is recommended to add the absolute path of the paired 
FASTQ files, but relative path should work also. Note, the delimiter is the tab (\t) character:
```
G15511	0	normal	normal_1	/samples/G15511/BL/C09DFACXX111207.1_1.fastq.gz	/samples/G15511/BL/C09DFACXX111207.1_2.fastq.gz
G15511	0	normal	normal_2	/samples/G15511/BL/C09DFACXX111207.2_1.fastq.gz	/samples/G15511/BL/C09DFACXX111207.2_2.fastq.gz
G15511	0	normal	normal_3	/samples/G15511/BL/C09DFACXX111207.3_1.fastq.gz	/samples/G15511/BL/C09DFACXX111207.3_2.fastq.gz
G15511	1	tumor	tumor_1	/samples/G15511/T/D0ENMACXX111207.1_1.fastq.gz	/samples/G15511/T/D0ENMACXX111207.1_2.fastq.gz
G15511	1	tumor	tumor_2	/samples/G15511/T/D0ENMACXX111207.2_1.fastq.gz	/samples/G15511/T/D0ENMACXX111207.2_2.fastq.gz
```

# Example TSV file for a normal/tumor pair with not-recalibrated BAM files
On the other hand, if you have pre-processed but not-recalibrated BAMs (that is the de-duplicated and realigned BAMs), their indexes and recalibration tables, you should use a structure like:
```
G15511	0	normal	pathToFiles/G15511.normal__1.md.real.bam	pathToFiles/G15511.normal__1.md.real.bai	pathToFiles/G15511.normal__1.recal.table
G15511	1	tumor	pathToFiles/G15511.tumor__1.md.real.bam	pathToFiles/G15511.tumor__1.md.real.bai	pathToFiles/G15511.tumor__1.recal.table
```
All the files will be created in the Preprocessing/NonRecalibrated/ directory, and by default a corresponding TSV file will also be deposited there. Generally, getting MuTect1 and Strelka calls on the preprocessed files should be done by:
```bash
nextflow run SciLifeLab/CAW --sample Preprocessing/NonRecalibrated/mysample.tsv --steps recalibrate,MuTect1,Strelka
```
# Example TSV file for a normal/tumor pair with recalibrated BAM files
The same way, if you have recalibrated BAMs (that is the de-duplicated and realigned and recalibrated BAMs) and their indexes, you should use a structure like:
```
G15511	0	normal	pathToFiles/G15511.normal__1.md.real.bam pathToFiles/G15511.normal__1.md.real.bai
G15511  1	tumor	pathToFiles/G15511.tumor__1.md.real.bam pathToFiles/G15511.tumor__1.md.real.bai
```
All the files will be in he Preprocessing/Recalibrated/ directory, and by default a corresponding TSV file will also be deposited there. Generally, getting MuTect1 and Strelka calls on the recalibrated files should be done by:

```bash
nextflow run SciLifeLab/CAW --sample Preprocessing/Recalibrated/mysample.tsv --steps skipPreprocessing,MuTect1,Strelka
```
