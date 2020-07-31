#! /usr/bin/env python

import sys, re, math, random

#VCF file whould be downloaded from ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp//release/20130502/ALL.wgs.phase3_shapeit2_mvncall_integrated_v5b.20130502.sites.vcf.gz

vcffile = "ALL.wgs.phase3_shapeit2_mvncall_integrated_v5b.20130502.sites.vcf"

outfile = "1000Genomes_20130502_SNP_maf0.3.vcf"
of=open(outfile, 'w')
for line in open(vcffile, 'r'):
    line=line.strip()
    if not line.startswith("#"):
        pma=re.compile('MULTI_ALLELIC')
        ma=pma.search(line)
        if not (ma):
            psnp=re.compile('VT=SNP')
            snp=psnp.search(line)
            if(snp):
                info=line.split("\t")[7]
                af_info=info.split(";")[1]
                fq=float(af_info.split("=")[1])
                if fq > 0.3:
                    of.write("%s\t%s\n" %(line.split("\t")[0], line.split("\t")[1]))

