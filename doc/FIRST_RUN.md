# Try the workflow
This small tutorial will explain to you how to run CAW on a small sample test data. Some variables are specific to Swedish UPPMAX cluster, but can be easily modified to suit any clusters.

## Make a test directory
```bash
mkdir test_CAW
cd test_CAW
```

## Copy and extract the sample test file
```bash
wget https://github.com/SciLifeLab/CAW/blob/master/data/tiny/tiny.tar.gz?raw=true -O tiny.tar.gz
tar -xvzf tiny.tar.gz
rm tiny.tar.gz
```

## Run the workflow
This workflow itself needs no installation. Nextflow will automatically fetch it from GitHub when launched if `SciLifeLab/CAW` is specified as the workflow name.
```bash
nextflow run SciLifeLab/CAW --sample tiny.tsv
```
If you're using a Swedish UPPMAX cluster, don't forget to provide your project ID.
```bash
nextflow run SciLifeLab/CAW --sample tiny.tsv --project <UPPMAX_project_ID>
```

# Other possibility for advance users

## Clone the repository and run the workflow
You can download the repository yourself from GitHub and run them directly:
```bash
git clone https://github.com/SciLifeLab/CAW
cd CAW
nextflow run main.nf --sample data/tsv/tiny.tsv --steps preprocessing
```

## Load Nextflow
If you're running on a Swedish UPPMAX cluster you can load Nextflow as an environment module:
```bash
module load Nextflow
```
Environnement variables are set up each time the module is loaded, so you might want to set them up after loading the module.
