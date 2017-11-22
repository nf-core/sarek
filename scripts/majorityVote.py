#!/usr/bin/env python

"""
By providing N distinct VCF files this utility generates N new VCFs containing calls by N, N-1, N-2... majority votes.
For example, if you have a set of calls like mutect.vcf strelka.vcf freebayes.vcf , it will give you three files like:
    callsN.vcf      - records where all the callers agree
    callsN-1.vcf    - records where all but one callers agree
    callsN-2.vcf    - at least one caller gives a call (union of records)
Usage:
    majorityVote.py -v set1.vcf,set2.vcf,set3.vcf,... -e etalon.vcf
    etalon.vcf contains the expected calls
"""

import sys
import vcf
import click

from sets import Set

@click.command(context_settings = dict( help_option_names = ['-h', '--help'] ))
@click.option('--vcf',  '-v', type=str, help='List of VCF files to merge', required=True)
@click.option('--etalon',  '-e', type=str, help='Etalon VCF to compare to ', required=False)
def mergeVCFs(vcf,etalon):
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
    # now time to print out the calls
    printCalls(votes,vcfSets,vcfFiles)
    # if there is an etalon, calculate concordance
    if etalon is not None:
        print "Comparing to etalon"
        etalonCalls = getVCFRecords(etalon)
        compareToEtalon(etalonCalls,votes)
    
def compareToEtalon(eCalls,votes):
    # go through each set, and compare to etalon calls
    for vSet in votes:
        concordanceCount = 0
        for r in vSet:
            if r in eCalls:
                concordanceCount = concordanceCount +1
        print "Concordance with votes " + str(votes.index(vSet)+1)+ " " + str(concordanceCount)


def printCalls(votes,vcfSets,vcfFiles):
    templateReader = vcf.Reader(open(vcfFiles[0], 'r'))

    # create writers, use the first input vcf as a template
    # number of output files:
    numberOfMajorityFiles = len(votes)
    mWriters = []
    for m in range(1,numberOfMajorityFiles+1):
        outFile = vcf.Writer(open("calls_" + str(m) +".vcf","w"), templateReader, lineterminator='\n')
        mWriters.append(outFile)
    
    # store all the records in a dictionary: of course there will be clashes,
    # but we want to go through the files only once
    allRecords = {}
    for f in vcfFiles:
        reader = vcf.Reader(open(f, 'r'))
        for record in reader:
            theKey = str(record.CHROM) + str(record.POS)
            allRecords[theKey] = record

    # now go through votes
    for sets in votes:
        idx = votes.index(sets)
        print "Writing sets with " + str(idx+1) + " votes to " + mWriters[idx].stream.name
        for v in sets:
            vl = list(v)
            voteKey = str(vl[0]) + str(vl[1])
            r = allRecords[voteKey]
            mWriters[idx].write_record(r)


def getVCFRecords(calls):
    """
    Store VCF record CROM,POS,ALT in a set
    """
    c = Set()
    reader = vcf.Reader(open(calls, 'r'))
    for record in reader:
        #c.add( (record.CHROM, record.POS, record.REF, str(record.ALT[0])) )
        c.add( (record.CHROM, record.POS) )
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
