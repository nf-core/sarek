#!/bin/bash
#SBATCH -A sens2016004
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
PREFIX=${2}_${DATE}
ln -fs /castor/project/proj_nobackup/CAW/containers containers

# nextflow specific stuff to save everything on /scratch

export NXF_TEMP=/scratch
export NXF_LAUNCHBASE=/scratch
export NXF_WORK=`pwd`/work
export NXF_HOME=/castor/project/proj_nobackup/nextflow
export PATH=${NXF_HOME}/bin:${PATH}
nextflow run ${CAW}/main.nf --step skipPreprocessing -profile singularityLocal --tools $2 --sample $1 -with-timeline ${PREFIX}.timeline.html -with-trace ${PREFIX}.trace.txt

# for annotation run
# sbatch -J sample-Ann ./whatever.sh result.vcf[.gz] snpEff
# nextflow run ${CAW}/main.nf --step annotate -profile singularityLocal --tools $2 --annotateVCF $1 --noReports --sample Preprocessing/Recalibrated/recalibrated.tsv -with-timeline ${PREFIX}.timeline.html -with-trace -with-trace ${PREFIX}.trace.txt
