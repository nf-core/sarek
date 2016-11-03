# Install and execute workflow
This small tutorial will explain to you how to install Nextflow and run CAW on a small sample test data.

Some variables are specific to Swedish UPPMAX cluster, but can be easily modified to suit any clusters.

## Install Nextflow
To use this pipeline, you need to have a working version of Nextflow installed. You can find more information about this pipeline tool at [nextflow.io](http://www.nextflow.io/). The typical installation of Nextflow looks like this:
```bash
curl -fsSL get.nextflow.io | bash
mv ./nextflow ~/bin
```
`~/bin` should be in your `$PATH`.

## Create Nextflow specific directories
The second one might have already been created when you installed Nextflow.
```bash
mkdir $HOME/glob/nxftmp
mkdir $HOME/.nextflow
```

## Configure environnement variables
Add to your `.bashrc`
```bash
export NXF_HOME=$HOME/.nextflow
export NXF_TEMP=${SNIC_TMP:-$HOME/glob/nxftmp}
export NXF_LAUNCHBASE=${SNIC_TMP:-$HOME/glob/nxftmp}
export NXF_WORK=$HOME/glob/work
export NXF_OPTS='-Xms1g -Xmx4g'
```

# Install and try the workflow

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
This workflow itself needs no installation - Nextflow will automatically fetch it from GitHub when run if `SciLifeLab/CAW` is specified as the workflow name.
```bash
nextflow run SciLifeLab/CAW --sample tiny.tsv --steps preprocessing --project <UPPMAX_project_ID>
```

# Other possibility for advance users

## Load Nextflow
If you're running on a Swedish UPPMAX cluster you can load Nextflow as an environment module:
```bash
module load Nextflow
```

## Workflow installation
You can download the repository yourself from GitHub and run them directly:
```bash
git clone https://github.com/SciLifeLab/CAW
cd CAW
nextflow run main.nf --sample data/tsv/tiny.tsv --steps preprocessing
```
