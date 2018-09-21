# How to run Sarek

This guide will take you through your first run of Sarek.
It is divided into two steps corresponding to the two main types of analysis offered by Sarek:
 - Run a Germline Analysis
 - Run a Somatic Analysis

This guide assumes you have internet access on the server where the analysis will take place. If you do not have that, please look into the [installation instructions](INSTALL_BIANCA.md) for the restricted access server Bianca at Uppmax, which should give an idea on how to adjust the following examples accordingly.


It is recommended to run Sarek within a [screen](https://www.gnu.org/software/screen/) or [tmux](https://tmux.github.io/) session.
This helps Sarek run uninterrupted until the analysis has finished.
Furthermore, Sarek is designed to be run on a single sample for a germline analysis or a set of samples from the same individual for a somatic analysis.
If more than one individual will be analysed, it is recommended that this is done in separate directories which is analysed separately.


## Update to latest version

To make sure that you have the latest version of Sarek, use:

```bash
nextflow pull SciLifeLab/Sarek
```

## Run the latest version

If there is a feature or bugfix you want to use in a resumed or re-analyzed run, you have to update the workflow to the latest version.
By default it is not updated automatically, so use something like:

```bash
nextflow run -latest SciLifeLab/Sarek/main.nf ... -resume
```

## Not on Uppmax
The commands used in this guide is suitable on how to run on a cluster at Uppmax.
To run these examples on a different infrastructure, there are a few things that needs to be changed.

 - Most likely, the `slurm` profile is not suitable to use.
 Find a more suitable one (or design your own) using the [configuration documentation](CONFIG.md)
 - The path for where reference genomes are located (specified in the `--genome_base` parameter) need to be modified.
 Use the instructions in the [reference documentation](REFERENCES.md) to make sure all the reference files are available.


## Run a Germline Analysis
This section presents a complete instruction to run a germline analysis using Sarek on a single sample.
Sarek will start the analysis by parsing a supplied input file in TSV format.
This file contains all the necessary information about the data and for the germline analysis it should have at least one line.
For more detailed information about how to construct TSV files for custom data, see [input documentation](INPUT.md).

For example, the file can be called `samples_germline.tsv` with the content (corresponding to columns: `subject gender status sample lane fastq1 fastq2`):

```
SUBJECT_ID  XX    0    SAMPLEID    1    /samples/normal_1.fastq.gz    /samples/normal_2.fastq.gz
```

The first workflow that will be run is contained in the `main.nf` file and performs the preprocessing step consisting of mapping, marking of duplicates and base recalibration. Running this command will launch a nextflow process in the terminal which in turn submits jobs (processes) to the SLURM queue.
```
nextflow run SciLifeLab/Sarek/main.nf \
--sample samples_germline.tsv \
-profile slurm  \
--project <your uppmax project id> \
--genome_base /sw/data/uppnex/ToolBox/hg38bundle \
--genome GRCh38
```

When the workflow has finished successfully it should print something similar to this:
```
Completed at: Fri Aug 31 05:10:07 CEST 2018
Duration    : 1d 13h 24m 51s
Success     : true
Exit status : 0
```
Make sure to check that the output states `Success : true` and not `Success : false`.
The results of the first step is located in the `Preprocessing` directory.
These files will be used in the next step, where the actual variant calling takes place.
Among other things, the preprocessing step should have created a new TSV file which is intended to be used as input for the variant calling step:
```
nextflow run SciLifeLab/Sarek/germlineVC.nf \
--sample Preprocessing/Recalibrated/recalibrated.tsv \
-profile slurm  \
--project <your uppmax project id> \
--genome_base /sw/data/uppnex/ToolBox/hg38bundle \
--genome GRCh38 \
--tools HaplotypeCaller
```
When successful (`Success : true`), this step should produce vcf file(s) within a `VariantCalling` directory.
The next workflow will annotate the found variants.
It is possible to specify the tools used for annotation (here VEP) and the variant-calling tools to use as input for annotation (here HaplotypeCaller).
```
nextflow run SciLifeLab/Sarek/annotate.nf \
--annotateTools HaplotypeCaller \
-profile slurm \
--project <your uppmax project id> \
--genome_base ~/Sarek/References/smallGRCh37 \
--tools VEP
```

Finally, run MultiQC to get an easily accessible report of all your analysis.
```
nextflow run SciLifeLab/Sarek/runMultiQC.nf \
-profile slurm
--project <your uppmax project id> \
```
## Run a Somatic Analysis

This section presents a complete instruction on how to run a somatic analysis using Sarek on two  samples from the same individual. In this case one normal sample and one tumour sample will be used. However, Sarek can also accept more than one tumour sample (i.e. relapses) for the same individual.

Note: Four out of five of the steps included in this example are identical or very similar to the steps included in the germline analysis example. Therefore, much of the information in this example is redundant compared to the first example.

Sarek will start the analysis by parsing a supplied input file in TSV format.
This file contains all the necessary information about the data and for the somatic analysis it should have at least two lines.
These lines have columns corresonding to `subject gender status sample lane fastq1 fastq2`.
For more detailed information about how to construct TSV files for custom data, see [input documentation](INPUT.md).

For example, the file can be called `samples_somatic.tsv` with the content:

```
SUBJECT_ID  XX    0    SAMPLEID1    1    /samples/normal_1.fastq.gz    /samples/normal_2.fastq.gz
SUBJECT_ID  XX    1    SAMPLEID2    1    /samples/tumour_1.fastq.gz    /samples/tumour_2.fastq.gz
```
The first workflow that will be run is contained in the `main.nf` file and performs the preprocessing step consisting of mapping, marking of duplicates and base recalibration. Running this command will launch a nextflow process in the terminal which in turn submits jobs (processes) to the SLURM queue.
```
nextflow run SciLifeLab/Sarek/main.nf \
--sample samples_somatic.tsv \
-profile slurm  \
--project <your uppmax project id> \
--genome_base /sw/data/uppnex/ToolBox/hg38bundle \
--genome GRCh38
```

When the workflow has finished successfully it should print something similar to this:
```
Completed at: Fri Aug 31 05:10:07 CEST 2018
Duration    : 1d 13h 24m 51s
Success     : true
Exit status : 0
```

Make sure to check that the output states `Success : true` and not `Success : false`.
The results of the first step is located in the `Preprocessing` directory.
These files will be used in the next two steps, where the actual variant calling takes place.
Among other things, the preprocessing step should have created a new TSV file which is intended to be used as input for the variant calling steps:

```
nextflow run SciLifeLab/Sarek/germlineVC.nf \
--sample Preprocessing/Recalibrated/recalibrated.tsv \
-profile slurm  \
--project <your uppmax project id> \
--genome_base /sw/data/uppnex/ToolBox/hg38bundle \
--genome GRCh38 \
--tools HaplotypeCaller
```
When successful (`Success : true`), this step should produce vcf file(s) within a `VariantCalling` directory.
The first variant calling step is actually the one from the germline analysis.
This is included here since information regarding germline variants is still useful for analysis of somatic variants.
The next variant calling step is the somatic specific analysis:
```
nextflow run SciLifeLab/Sarek/somaticVC.nf \
--sample Preprocessing/Recalibrated/recalibrated.tsv \
-profile slurm \
--project <your uppmax project id> \
--genome_base /sw/data/uppnex/ToolBox/hg38bundle \
--genome GRCh38 \
--tools Strelka
```
When successful (`Success : true`), this step should produce vcf file(s) within the `VariantCalling` directory separate from the germline vcf file.
The next workflow will annotate the found variants.
It is possible to specify the tools used for annotation (here VEP) and the variant-calling tools to use as input for annotation (here HaplotypeCaller and Strelka).
```
nextflow run SciLifeLab/Sarek/annotate.nf \
--annotateTools HaplotypeCaller,Strelka \
-profile slurm \
--project <your uppmax project id> \
--genome_base ~/Sarek/References/smallGRCh37 \
--containerPath \
--tools VEP
```

Finally, run MultiQC to get an easily accessible report of all your analysis.
```
nextflow run SciLifeLab/Sarek/runMultiQC.nf \
-profile slurm
--project <your uppmax project id> \
```
