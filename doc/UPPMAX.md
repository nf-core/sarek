# Install and execute workflow
This small tutorial will explain to you how to install Nextflow and run CAW on a small sample test data.

Some variables are specific to Swedish UPPMAX cluster, but can be easily modified to suit any clusters.

## Install Nextflow
To use this pipeline, you need to have a working version of Nextflow installed. You can find more information about this pipeline tool at [nextflow.io](http://www.nextflow.io/). The typical installation of Nextflow looks like this:
```bash
curl -fsSL get.nextflow.io | bash
mv ./nextflow ~/bin
```

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

## Download the workflow config file
The config file is based on a [config file](https://raw.githubusercontent.com/SciLifeLab/CAW/master/config/milou.config) specific to Swedish UPPMAX milou cluster, but can be easily modified to suit any clusters.
```bash
wget https://raw.githubusercontent.com/SciLifeLab/CAW/master/config/milou.config -O $NXF_HOME/config
```
If you're using this config file, don't forget to edit the line `'-A b2015110'` to contain your own UPPMAX project identifier instead.

## Copy the sample test file
```bash
wget https://raw.githubusercontent.com/SciLifeLab/CAW/master/data/tsv/tiny-github.tsv
```

## Run the workflow
This workflow itself needs no installation - Nextflow will automatically fetch it from GitHub when run if `SciLifeLab/CAW` is specified as the workflow name.
```bash
nextflow run SciLifeLab/CAW --sample tiny-github.tsv --steps preprocessing
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
nextflow run CAW/main.nf
```

It's possible to run the full workflow on a data set specified in the tsv file. For example:
```bash
nextflow run SciLifeLab/CAW -c config/milou.config --sample data/tsv/tiny.tsv --steps preprocessing
```
will run the workflow on a small testing dataset. To try on your own data, you just have to edit your own tsv file.
