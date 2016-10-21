# UPPMAX

## Loading Nextflow
If you're running on a Swedish UPPMAX cluster you can load Nextflow as an environment module:
```bash
module load Nextflow
```

## Configure Nextflow
I would advice you to create a Nextflow directory in your home directory and use it as `NXF_HOME`. That way, when launching a Nextflow script directly from github, Nextflow will be able to download it.
```bash
mkdir ~/.nextflow
```
And then adding the following line to your `~/.bashrc` or `~/.bash_profile` file:
```bash
export NXF_HOME=~/.nextflow
```
Every time you load the Nextflow module, you will get a warning about setting environment variables. Adding the following lines to your your `~/.bashrc` or `~/.bash_profile` file, will set these automatically at login, and should get rid of the warning:
```bash
export NXF_LAUNCHBASE=$SNIC_TMP
export NXF_TEMP=$SNIC_TMP
```
To make sure Nextflow doesn't eat all memory, you an also add the following line to your `~/.bashrc` or `~/.bash_profile` file:
```bash
export NXF_OPTS='-Xms1g -Xmx4g '
```

## Configure the workflow for easy execution
Next, you will need to set up a config file so that Nextflow knows how to run the workflow and paths to reference files. You can find an example configuration file for UPPMAX (milou) with this repository:
[`milou-github.config`](https://github.com/SciLifeLab/CAW/blob/master/config/milou-github.config).

Copy this file in your `$NXF_HOME` which should be `~/.nextflow` af `config` and edit the line `'-A b2015110'` to contain your own UPPMAX project identifier instead.

```bash
wget https://raw.githubusercontent.com/SciLifeLab/CAW/master/config/milou-github.config $NXF_HOME/config
```
