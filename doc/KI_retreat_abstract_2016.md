### Cancer Analysis Workflow Of Tumor/Normal Pairs At The National Genomics Infrastructure Of SciLifeLab

One of the most prominent usage of NGS is whole genome sequencing (WGS). At
SciLifeLab with the available infrastructure it is possible to produce data for
normal and tumor sample pairs in a couple of days. However, building a robust
and reliable bioinformatics workflow to find somatic mutations is challenging:
tumor samples are heterogeneous, likely contain structural variants and
multiple sub-clones besides the normal tissue.

We are presenting our workflow that is designed to analyze WGS tumor/normal
data in a high-throughput environment. The framework is based on the nextflow
domain-specific language on top of java/groovy. Using nextflow we are able to
utilize both the slurm load balancing environment and local execution,
implementing data flow forks and joins, call external software etc. Individual
sub-steps of a complex flow can be connected and restarted after failure from
the last execution point.

The actual preprocessing workflow is based on BWA as an aligner and GATK best
practice steps. To achieve a consensus variant call different variant callers
can be added, currently MuTect2, Strelka and VarDict are supported, more to be
added. Structural variants are going to be estimated by manta, ploidy and
sample heterogeneity is measured by ASCAT. The expected output of the workflow
is a VCF file presenting filtered, prioritized and annotated polymorphisms.

As the nextflow environment is flexible, we can add other tools or remove
obsolete ones as development progresses. The goal is to build a workflow for
cancer genome analysis that can be deployed to both research and clinical
environments and are going to be included as a standard workflow at NGI during
the fall of 2016.

