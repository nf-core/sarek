#!/bin/bash
#SBATCH -p node
#SBATCH -t 168:00:00

set -xeuo pipefail

# skeleton script to launch CAW jobs on scratch with sbatch on bianca

GENOME=GRCh38
SAMPLE=''
STEP=''
TOOLS=false

while [[ $# -gt 0 ]]
do
  key=$1
  case $key in
    -g|--genome)
    GENOME=$2
    shift # past argument
    shift # past value
    ;;
    -s|--sample)
    SAMPLE=$2
    shift # past argument
    shift # past value
    ;;
    --step)
    STEP=$2
    shift # past argument
    shift # past value
    ;;
    -t|--tools)
    TOOLS=$2
    shift # past argument
    shift # past value
    ;;
    *) # unknown option
    shift # past argument
    ;;
  esac
done

# Make a specific prefix for logs (trace + timeline)
DATE=`date +%Y-%b-%d-%H%M`
if [[ $TOOLS ]]
then
  PREFIX=${TOOLS}_${DATE}
else
  PREFIX=${DATE}
fi

# Use the default deployed CAW version.
# Other versions (including testing) are at /castor/project/proj_nobackup/CAW/
CAW=/castor/project/proj_nobackup/CAW/default

# Link containers
ln -fs /castor/project/proj_nobackup/CAW/containers containers

# Configure Nextflow to save everything on /scratch
export NXF_TEMP=/scratch
export NXF_LAUNCHBASE=/scratch
export NXF_WORK=/scratch
export NXF_HOME=/castor/project/proj_nobackup/nextflow
export PATH=${NXF_HOME}/bin:${PATH}

# save as whatever, and launch like:
# sbatch -A [project name] -J sample-MuTect2 ./whatever.sh sample.tsv MuTect2

function run_caw() {
  nextflow run ${CAW}/main.nf $@ -with-timeline ${PREFIX}.timeline.html -with-trace ${PREFIX}.trace.txt
}

# MAPPING
if [[ $STEP == MAPPING ]]
then
  run_caw() --sample ${SAMPLE} --step mapping
fi

# RUN A SPECIFIC TOOL:
#  -s or --sample the sample TSV file
#  -t or --tools the tool to run
# Don't forget to specify the right SNIC project
# or add as the second line:
#SBATCH -A [project name]

if [[ $STEP == VARIANTCALLING ]]
then
  run_caw() --sample ${SAMPLE} --step variantcalling --tools ${TOOLS}
fi

# for annotation run
# sbatch -A [project name] -J sample-Ann ./whatever.sh result.vcf[.gz] snpEff
# nextflow run ${CAW}/main.nf --step annotate --tools $2 --annotateVCF $1 --noReports --sample Preprocessing/Recalibrated/recalibrated.tsv -with-timeline ${PREFIX}.timeline.html -with-trace -with-trace ${PREFIX}.trace.txt

if [[ $STEP == ANNOTATE ]]
then
  run_caw() --step annotate --tools ${TOOLS} --annotateVCF ${SAMPLE}
fi
