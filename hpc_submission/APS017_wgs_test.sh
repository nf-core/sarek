#!/bin/bash
#SBATCH --partition=componc_cpu,componc_gpu
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --time=48:00:00
#SBATCH --mem=8GB
#SBATCH --job-name=sarek
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=preskaa@mskcc.org
#SBATCH --output=slurm%j_TCDO-SAR-061.out


## activate nf-core conda environment
source /home/preskaa/miniconda3/bin/activate nf-core

## load modules
module load singularity/3.7.1
module load java/20.0.1
## example samplesheet
## technical replicates get merged ...
samplesheet=${HOME}/sarek/resources/wgs_test_samplesheet.csv
## specify path to out directory
outdir=/data1/shahs3/users/preskaa/APS017_Archive/sarek

## reference genome for bwacmp (if sample is PDX)
mouse_refgenome=/data1/shahs3/isabl_data_lake/assemblies/WGS-MM10/mouse/mm10_build38_mouse.fasta
### human reference genome
refgenome=/data1/shahs3/isabl_data_lake/assemblies/GRCh38-P14/GRCh38.primary_assembly.genome.fa
ref_index=/data1/shahs3/isabl_data_lake/assemblies/GRCh38-P14/GRCh38.primary_assembly.genome.fa.fai

## last two flags trigger chopper to differentiate mouse from human reads for PDX samples
## these flags should not be used for human samples
cd ${outdir}

nextflow run apsteinberg/sarek \
  -c ${HOME}/nanoseq/conf/iris.config \
  -profile singularity,slurm \
  --input ${samplesheet} \
  --outdir ${outdir} \
  -work-dir ${outdir}/work \
  --trim_fastq \
  --aligner bwa-mem \
  --save_mapped \
  --save_output_as_bam \
  --fasta ${refgenome} \
  --fasta_fai ${ref_index} \
  --email preskaa@mskcc.org


#nextflow run apsteinberg/nanoseq -resume 6c03bf60-99ea-41cd-a949-c30986899f14
