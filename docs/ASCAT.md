# ASCAT

## Introduction

Ascat is a software for performing allele-specific copy number analysis of tumor samples and for estimating tumor ploidy and purity (normal contamination). Ascat is written in R and available here: [github.com/Crick-CancerGenomics/ascat](https://github.com/Crick-CancerGenomics/ascat).

To run Ascat on NGS data we need .bam files for the tumor and normal samples, as well as a loci file with SNP positions. If Ascat is run on SNP array data, the loci file contains the SNPs on the chip. When runnig Ascat on NGS data we can use the same loci file, for exampe the one corresponding to the AffymetrixGenome-Wide Human SNP Array 6.0, but we can also choose a loci file of our choice with i.e. SNPs detected in the 1000 Genomes project.

### BAF and LogR values

Running Ascat on NGS data requires that the .bam files are converted into BAF and LogR values. This can be done using the software [AlleleCount](https://github.com/cancerit/alleleCount) followed by a simple R script. AlleleCount extracts the number of reads in a bam file supporting each allele at specified SNP positions. Based on this, BAF and logR can be calculated for every SNP position i as:

```
BAFi(tumor)=countsBi(tumor)/(countsAi(tumor)+countsBi(tumor))
BAFi(normal)=countsBi(normal)/(countsAi(normal)+countsBi(normal))
LogRi(tumor)=log2((countsAi(tumor)+countsBi(tumor))/(countsAi(normal)+countsBi(normal)) - median(log2((countsA(tumor)+countsB(tumor))/(countsA(normal)+countsB(normal)))
LogRi(normal)=0
```

For male samples, the X chromosome markers have special treatment:

```
LogRi(tumor)=log2((countsAi(tumor)+countsBi(tumor))/(countsAi(normal)+countsBi(normal))-1 - median(log2((countsA(tumor)+countsB(tumor))/(countsA(normal)+countsB(normal))-1)
```

where:

- `i` corresponds to the postions of all SNPs in the loci file.
- `CountsA` and `CountsB` are vectors containing number of reads supporting the `A` and `B` alleles of all SNPs
- `A` = the major allele
- `B` = the minor allele
- Minor and major alleles are defined in the loci file (it actually doesn't matter which one is defied as A and B in this application).

Calculation of LogR and BAF based on AlleleCount output is done as in [runASCAT.R](https://github.com/cancerit/ascatNgs/tree/dev/perl/share/ascat/runASCAT.R)

### Loci file

The loci file was created based on the 1000Genomes latest release (phase 3, releasedate 20130502), available [here](ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp//release/20130502/ALL.wgs.phase3_shapeit2_mvncall_integrated_v5b.20130502.sites.vcf.gz). The following filter was applied: Only bi-allelc SNPs with minor allele frequencies > 0.3. The filtered file can be found on [export.uppmax.uu.se](https://export.uppmax.uu.se/b2015110/caw-references/b37/1000G_phase3_20130502_SNP_maf0.3.loci.tar.bz2) and is stored on Milou in:

```
/sw/data/uppnex/ToolBox/ReferenceAssemblies/hg38make/bundle/2.8/b37/1000G_phase3_20130502_SNP_maf0.3.loci
```

The loci file was originally generated for GRCh37. It was translated into GRCh38 using the tool liftOver available at the UCSC Genome Browser. To run liftOver the loci file was first written in bed format:

```
awk '{print "chr"$1":"$2"-"$2}' 1000G_phase3_20130502_SNP_maf0.3.loci > 1000G_phase3_20130502_SNP_maf0.3.bed
```

Using the web interface to liftOver at [genome.ucsc.edu](https://genome.ucsc.edu/cgi-bin/hgLiftOver) the file was translated into GRCh38 coordinates. LiftOver was possible for 3261270 out of 3268043 SNPs. The converted SNP positions were printed in the format required by AlleleCounter by:

```
more hglft_genome_5834_13aba0.bed | awk 'BEGIN{FS="chr"} {print $2}' | awk 'BEGIN{FS="-"} {print $1}' | awk 'BEGIN{FS=":";OFS="\t"} {print $1,$2}' > 1000G_phase3_GRCh38_maf0.3.loci
```

The loci file in GRCh38 coordinates is stored on Milou in:

```
/sw/data/uppnex/ToolBox/ReferenceAssemblies/hg38make/bundle/2.8/1000G_phase3_GRCh38_maf0.3.loci
```

## Running on Milou

### Run AlleleCount

AlleleCount is installed as part of the `bioinfo-tools` module on Milou. It runs on single bam files (tumor and normal separately) with the command below:

```bash
$ module load bioinfo-tools alleleCount
$ alleleCounter -l /sw/data/uppnex/ToolBox/ReferenceAssemblies/hg38make/bundle/2.8/b37/1000G_phase3_20130502_SNP_maf0.3.loci -r /sw/data/uppnex/ToolBox/ReferenceAssemblies/hg38make/bundle/2.8/b37/human_g1k_v37_decoy.fasta -b sample.bam -o sample.allecount
```

### Convert allele counts to LogR and BAF values

The allele counts can then be converted into LogR and BAF values using the script `convertAlleleCounts.r`. Usage for a male sample (`Gender = "XY"`, replace with `Gender = "XX"` for a female sample):

```bash
sbatch -A PROJID -p core -n 1 -t 240:00:00 -J convertAllelecounts -e convertAllelecounts.err -o convertAllelecounts.out /path/to/your/Sarek-fork/convertAlleleCounts.r tumor_sample tumor.allelecount normal_sample normal.allelecount XY
```

This creates the BAF and LogR data for the tumor and normal samples, to be used as input to ASCAT.

### Run ASCAT

The script "run_ascat.r" can be used to run ASCAT in the simplest possible way without compensating for the local CG content across the genome. It calls the main ASCAT R script [ascat.R](https://github.com/Crick-CancerGenomics/ascat/tree/master/ASCAT/R/ascat.R).

```bash
sbatch -A PROJID -p core -n 1 -t 240:00:00 -J ascat -e ascat.err -o ascat.out run_ascat.r tumor_baf tumor_logr normal_baf normal_logr
```

## Flowchart

![Overview of ASCAT process](images/ascat.jpg "ASCAT")

--------------------------------------------------------------------------------

[![](images/SciLifeLab_logo.png "SciLifeLab")][scilifelab-link]
[![](images/NGI_logo.png "NGI")][ngi-link]
[![](images/NBIS_logo.png "NBIS")][nbis-link]

[nbis-link]: https://www.nbis.se/
[ngi-link]: https://ngisweden.scilifelab.se/
[scilifelab-link]: https://www.scilifelab.se/
