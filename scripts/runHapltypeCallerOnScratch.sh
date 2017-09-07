#!/bin/bash -xl
#SBATCH -A b2011185
#SBATCH -p node
#SBATCH -t 168:00:00
module load Nextflow
export NXF_TEMP=/scratch
export NXF_LAUNCHBASE=/scratch
export NXF_WORK=/scratch
nextflow run /proj/b2011185/nobackup/private/analysis/CAW/CAW/main.nf --project b2011185 --step skippreprocessing --tools haplotypecaller --genome GRCh37 --sample $1
