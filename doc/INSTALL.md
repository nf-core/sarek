# How to install the workflow from scratch on a Linux machine

## Install Lmod

To control which version is used for each tools and dependencies for the workflow, a module environnement system as lmod is used.

Follow instructions for installing it on [lmod.readthedocs.io](http://lmod.readthedocs.io/en/latest/index.html)

## Install slurm

The idea here is to emulate a slurm ordonnancing system on a non-cluster machine to run the workflow.

Enable
OPTIONS="--force --key-file /etc/munge/munged.key --num-threads 1"
on '/etc/default/munge'

## Tools and dependencies
https://sourceforge.net/projects/bio-bwa/files/
- bwa 0.7.8
https://software.broadinstitute.org/gatk/download/index
- GATK 3.3-0
- GATK 3.6 [MuTect2]
https://sites.google.com/site/strelkasomaticvariantcaller/home/download
- Strelka 1.0.15
- gcc 4.9.2
- java sun_jdk 1.8.0_92
- manta 0.27.1
- MuTect1 1.1.5
- Nextflow 0.17.3
- perl 5.18.4
- picard 1.118
- R 3.2.3
- samtools 0.1.19 [Manta]
- samtools 1.3
- snpeff 4.2
- vardict 1.4.5