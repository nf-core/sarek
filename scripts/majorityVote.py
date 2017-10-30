#!/usr/bin/env python

"""
By providing N distinct VCF files this utility generates N new VCFs containing calls by N, N-1, N-2... majority votes.
For example, if you have a set of calls like mutect.vcf strelka.vcf freebayes.vcf , it will give you three files like:
    callsN.vcf      - records where all the callers agree
    callsN-1.vcf    - records where all but one callers agree
    callsN-2.vcf    - at least one caller gives a call (union of records)
Usage:
    majorityVote.py -v set1.vcf,set2.vcf,set3.vcf,...
"""


import vcf
import click

from sets import Set

@click.command(context_settings = dict( help_option_names = ['-h', '--help'] ))
@click.option('--vcf',  '-v', type=str, help='List of VCF files to merge', required=True)
def mergeVCFs(vcf):
    vcfFiles = vcf.split(",")
    vcfSets = []
    print "Processing VCF files: "
    for vcfCalls in vcfFiles:
        print "\t" + vcfCalls
        vcfSets.append( getVCFRecords(vcfCalls) )
    # OK, now print out the results
    # one day I will write a proper VCF exporter, but not today
    print "#callers ##calls"
    votes = calculateVotes(vcfSets)
    n = 1
    for v in votes:
        print n,len(v)
        n = n+1

def getVCFRecords(calls):
    """
    Store VCF record CROM,POS,ALT in a set
    """
    c = Set()
    reader = vcf.Reader(open(calls, 'r'))
    for record in reader:
        c.add( (record.CHROM, record.POS, record.REF, str(record.ALT[0])) )
    return c

def calculateVotes(sets):
    """
    Generates intersections and votes from call sets. 
    The list member votes[N-1] will contain the set of calls supported by at least N callers
    """
    votes = []
    for i in range(0,len(sets)):
        votes.append( atLeastN(i+1,sets) )
#    votes.reverse()
    return votes

def atLeastN(n,sets):
    # I know we are generating a union every time, but right now I do not want to optimize yet
    union = getUnion(sets)
    #print "Votes for " + str(n)
    callsInAtN = Set()          # this will contain calls that are in called at least N callers
    for record in union:        # for every record in the union set
        count = 0               # count of votes
        for s in sets:          # pick a set
            if record in s:     # if the record (that is in the union) is also in the set
                count = count+1 # increment count
        # now we have the vote count for all call sets
        if count >= n:
            callsInAtN.add(record)
    #print callsInAtN
    return callsInAtN

def getUnion(sets):
    union = Set()
    for s in sets:
        union |= s
    return union

if __name__ == "__main__":
        mergeVCFs()
