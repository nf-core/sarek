#!/bin/bash -xl
#SBATCH -A sens2016004
#SBATCH -p node
#SBATCH -t 168:00:00
module load Nextflow
export NXF_TEMP=/scratch
export NXF_LAUNCHBASE=/scratch
export NXF_WORK=/scratch
nextflow run ${HOME}/CAW/main.nf -c ${HOME}/CAW/nextflow.config -profile localhost --project sens2016004 --step preprocessing --sample $1
