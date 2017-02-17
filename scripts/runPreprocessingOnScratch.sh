#!/bin/bash -xl
#SBATCH -A sens2016004
#SBATCH -p node
#SBATCH -t 168:00:00
module load Nextflow
export NXF_TEMP=/scratch
export NXF_LAUNCHBASE=/scratch
export NXF_WORK=/scratch
nextflow run /home/szilva/CAW/main.nf -c /home/szilva/CAW/nextflow.config -profile localhost --project sens2016004 --verbose --steps preprocessing --sample $1
