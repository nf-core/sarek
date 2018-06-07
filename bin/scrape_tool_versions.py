#!/usr/bin/env python
from __future__ import print_function
from collections import OrderedDict
import re

regexes = {
    'AlleleCount': ['v_allelecount.txt', r"(\S+)"],
    'ASCAT': ['v_ascat.txt', r"(\d\.\d+)"],
    'bcftools': ['v_bcftools.txt', r"bcftools (\S+)"],
    'BWA': ['v_bwa.txt', r"Version: (\S+)"],
    'FastQC': ['v_fastqc.txt', r"FastQC v(\S+)"],
    'GATK': ['v_gatk.txt', r"GATK version(\S+)"],
    'htslib': ['v_samtools.txt', r"htslib (\S+)"],
    'Manta': ['v_manta.txt', r"([0-9.]+)"],
    'MultiQC': ['v_multiqc.txt', r"multiqc, version (\S+)"],
    'Nextflow': ['v_nextflow.txt', r"(\S+)"],
    'FreeBayes': ['v_freebayes.txt', r"version:  v(\d\.\d\.\d+)"],
    'Picard': ['v_picard.txt', r"Picard version:(\d\.\d\.\d+)"],
    'Qualimap': ['v_qualimap.txt', r"QualiMap v.(\S+)"],
    'R': ['v_r.txt', r"R version (\S+)"],
    'samtools': ['v_samtools.txt', r"samtools (\S+)"],
    'Sarek': ['v_sarek.txt', r"(\S+)"],
    'SnpEff': ['v_snpeff.txt', r"version SnpEff (\S+)"],
    'Strelka': ['v_strelka.txt', r"([0-9.]+)"],
    'vcftools': ['v_vcftools.txt', r"([0-9.]+)"],
    'VEP': ['v_vep.txt', r"ensembl-vep          : (\S+)"],
}
results = OrderedDict()
results['Sarek'] = '<span style="color:#999999;\">N/A</span>'
results['Nextflow'] = '<span style="color:#999999;\">N/A</span>'
results['BWA'] = '<span style="color:#999999;\">N/A</span>'
results['samtools'] = '<span style="color:#999999;\">N/A</span>'
results['htslib'] = '<span style="color:#999999;\">N/A</span>'
results['GATK'] = '<span style="color:#999999;\">N/A</span>'
results['Picard'] = '<span style="color:#999999;\">N/A</span>'
results['Manta'] = '<span style="color:#999999;\">N/A</span>'
results['Strelka'] = '<span style="color:#999999;\">N/A</span>'
results['FreeBayes'] = '<span style="color:#999999;\">N/A</span>'
results['AlleleCount'] = '<span style="color:#999999;\">N/A</span>'
results['R'] = '<span style="color:#999999;\">N/A</span>'
results['ASCAT'] = '<span style="color:#999999;\">N/A</span>'
results['SnpEff'] = '<span style="color:#999999;\">N/A</span>'
results['VEP'] = '<span style="color:#999999;\">N/A</span>'
results['FastQC'] = '<span style="color:#999999;\">N/A</span>'
results['Qualimap'] = '<span style="color:#999999;\">N/A</span>'
results['bcftools'] = '<span style="color:#999999;\">N/A</span>'
results['vcftools'] = '<span style="color:#999999;\">N/A</span>'
results['MultiQC'] = '<span style="color:#999999;\">N/A</span>'

# Search each file using its regex
for k, v in regexes.items():
    try:
      with open(v[0]) as x:
        versions = x.read()
        match = re.search(v[1], versions)
      if match:
        results[k] = "v {}".format(match.group(1))
    except Exception as FileNotFoundError:
      print("No such file:", v[0])

# Remove empty keys (defining them above ensures correct order)
for k in ['Sarek', 'Nextflow', 'BWA', 'samtools', 'htslib', 'GATK', 'Picard', 'Manta', 'Strelka', 'FreeBayes', 'AlleleCount', 'R', 'ASCAT', 'SnpEff', 'VEP', 'FastQC', 'Qualimap', 'bcftools', 'vcftools', 'MultiQC']:
    if results[k] == '<span style="color:#999999;\">N/A</span>':
        del(results[k])

# Dump to YAML
print ('''
id: 'Sarek'
order: -1000
section_href: 'https://github.com/SciLifeLab/Sarek'
plot_type: 'html'
description: 'tool versions are collected at run time from output.'
data: |
  <dl class="dl-horizontal" style="margin-bottom:0;">
''')
for k,v in results.items():
    print("        <dt>{}</dt><dd>{}</dd>".format(k,v))
print ("    </dl>")
