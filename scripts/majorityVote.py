#!/usr/bin/env python
import vcf
import click
from vcf import utils
from sets import Set
"""
By providing N distinct VCF files this utility generates N new VCFs containing calls by N, N-1, N-2... majority votes.
For example, if you have a set of calls like mutect.vcf strelka.vcf freebayes.vcf , it will give you three files like:
    callsN.vcf      - records where all the callers agree
    callsN-1.vcf    - records where all but one callers agree
    callsN-2.vcf    - at least one caller gives a call (union of records)
The template VCF files is for a header
Usage:
    majorityVote.py -v set1.vcf,set2.vcf,set3.vcf,... -t template
"""


@click.command(context_settings = dict( help_option_names = ['-h', '--help'] ))
@click.option('--vcfRecords',  '-v', type=str, help='List of VCF files to merge', required=True)
@click.option('--template',  '-t', type=str, help='Template VCF with metadata (records ignored)', required=True)
def mergeVCFs(vcfrecords,template):
    vcfFiles = vcfrecords.split(",")
    # having N input VCF files we will need N writers as well
    writers = []
    for voteN in range(1,len(vcfFiles)+1):
        writers.append( vcf.Writer( open("calls_"+str(voteN)+".vcf","w"),vcf.Reader(open(template, 'r'))))
    # get readers for VCF files
    readers = []
    for r in vcfFiles:
        readers.append( vcf.Reader(open(r, 'r')))

    for records in utils.walk_together( *readers ):
        # count the number of non-empty records
        count = len(filter(None, records))
        toWrite = next(item for item in records if item is not None)
        for i in range(0,count):
            writers[i].write_record(toWrite)


if __name__ == "__main__":
        mergeVCFs()
