# Run Sarek using Conda

> /!\\
This functionality is still under development, so we would recommend to use it with care.

> /!\\
A conda environment is defined but does not contain all the tools needed for fully using Sarek.

To use Sarek with conda, first make sure that you have conda installed.
We recommend miniconda:
https://conda.io/miniconda.html

Sarek comes with a conda environment definition - a file called
[`environment.yml`](../environment.yml) which lists conda channels and package names/versions.

To run Sarek with Conda, you can use/modify the profile `conda`, or follow the following steps:

1. Create a new environment using the [`environment.yml`](../environment.yml) file:

```bash
# Download the environment.yml file
curl https://raw.githubusercontent.com/SciLifeLab/Sarek/master/environment.yml -o environment.yml

# Create a new conda environment using it
conda env create -f environment.yml
```

2. Create a nextflow config file as indicated below (Usually as `~/.nextflow/config` file):

```groovy
process {
  beforeScript = { 'module load anaconda; set +u; source activate sarek-2.2.2; set -u;' }
}
```

3. While launcing Sarek specify the nextflow config file crreated in step 2 using `-c` option:

```bash
nextflow run SciLifeLab/Sarek/main.nf \
  --sample samples_germline.tsv \
  -profile conda \
  --genome_base /sw/data/uppnex/ToolBox/hg38bundle \
  --genome GRCh38 \
  -c ~/.nextflow/config
```
