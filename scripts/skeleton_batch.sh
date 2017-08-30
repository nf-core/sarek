#!/bin/bash
#SBATCH -A sens2017102
#SBATCH -p node
#SBATCH -t 168:00:00
set -xeuo pipefail
# skeleton script to launch nextflow/singularity jobs with slurm on bianca
# save as whatever, and launch like:
# sbatch -J sample-MuTect2 ./whatever.sh sample.tsv MuTect2
# 
# PARAMETERS:
#  $1 is the sample TSV file
#  $2 the tool to run
#
# use the default deployed CAW version. "testing" versions should be at /castor/project/proj_nobackup/CAW/testing
CAW=/castor/project/proj_nobackup/CAW/default
DATE=`date +%Y-%b-%d-%H%M`
PREFIX=${1%.tsv}_${DATE}
TRACE=${PREFIX}.trace
TIMELINE=${PREFIX}.timeline.html

# nextflow specific stuff to save everything on /scratch

export NXF_TEMP=/scratch
export NXF_LAUNCHBASE=/scratch
export NXF_WORK=/scratch

nextflow run ${CAW}/main.nf -c ${CAW}/nextflow.config --step skipPreprocessing -with-singularity --profile singularityLocal --tools $2 $1
