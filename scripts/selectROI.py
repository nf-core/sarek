#!/usr/bin/env python
import os
import argparse
import subprocess

def parse_args():
    parser = argparse.ArgumentParser(description='Selecting Region Of Interest from BAM files')
    parser.add_argument('-v', help='List of VCF files to process like -v tumour.vcf,normal.vcf',required=False, dest="vcf")
    parser.add_argument('-a', help='The BAM files (indexed) to read from like -b tumour.bam,normal.bam', required=True, dest="bam")
    parser.add_argument('-b', help='The BED file with gene coordinates', required=False, dest="bed")
    parser.add_argument('-p', help='The prefix of the target BAM file to write to (postfix is the input filename)', required=False, dest="prefix", default="ROI.")

    return parser.parse_args()


def selectROI(vcf,bam,bed,prefix):
    # this is the set where loci are stored
    # VCF records are (chromosome,location) tuples
    # BEDs are (chromosome,start,end) tuples
    callDict = dict()

    print("Processing VCF files: ")
    vcfFiles = vcf.split(",")
    for vcfCalls in vcfFiles:
        print("\t" + vcfCalls)
        addVCFRecords(vcfCalls,callDict)
    print("Collected " + str(len(callDict)) + " chromosomes from VCFs")

    if bed:
        print("Processing BED files: ")
        bedFiles = bed.split(",")
        for bedLoci in bedFiles:
            print("\t" + bedLoci)
            addBEDRecords(bedLoci,callDict)
        print("Collected " + str(len(callDict)) + " records from BEDs")

    bamFiles = bam.split(",")
    print("Saving reads from BAMs:")
    saveCallsToBAM(callDict,bamFiles,prefix)

def saveCallsToBAM(callDict, bamFiles, prefix):
    for bam in bamFiles:
        print("\t" + bam)
        # getting reads one-by-one takes aages. We are going to try to use htslib/samtools
        # by making a command line to get ROI for each chromosome, and 
        # merge/sort the BAMs later
        filesToMerge = ""
        for chrom in callDict:
            cmdline = ""
            for gcoord in callDict[chrom]:
                if type(gcoord) is tuple:       # BED record in tuple
                    start = gcoord[0]
                    end = gcoord[1]
                else:                   # VCF record
                    start = gcoord - 500
                    end = gcoord - 500
                cmdline = cmdline + chrom + ":" + str(start) + "-" + str(end) + " "
            # let's write the file for the chromosome
            cmdline = cmdline.rstrip()
            chunkBAM = chrom + "_tmpROI.bam "
            filesToMerge = filesToMerge + chunkBAM
            subprocess.call("samtools view " + bam + " -bo " + chunkBAM + cmdline, shell=True)
            print("Wrote chunk for " + chrom + " to " + chunkBAM)
        subprocess.call("samtools merge -f merged.roi.bam " + filesToMerge, shell=True)
        target = os.path.dirname(bam) + "/" + prefix + "." + os.path.basename(bam)
        subprocess.call("samtools sort -@6 -o "+ target + " merged.roi.bam ", shell=True)
        print("Sorted file written to " + target)

def addVCFRecords(VCFFile,callDict):
    """
    Store VCF record CHROM,POS in a dict 
    """
    VCFlines = [line.rstrip() for line in open(VCFFile,'rt')]
    for variant in VCFlines:
        if not variant.startswith("#"):      # ignore header
            vcfCols = variant.split()
            if vcfCols[0] not in callDict:   # if chromosome not in dictionary yet
                callDict[vcfCols[0]] = []
            callDict[vcfCols[0]].append(int(vcfCols[1]))

def addBEDRecords(BEDFile,callDict):
    """
    Store BED record CHROM,START,END in a set
    """
    BEDlines = [line.rstrip() for line in open(BEDFile,'rt')]
    for locus in BEDlines:
        BEDrecord = locus.split()
        if BEDrecord[0] not in callDict:
            callDict[BEDrecord[0]] = []
        callDict[BEDrecord[0]].append( (int(BEDrecord[1]),int(BEDrecord[2])) )

if __name__ == "__main__":
    args = parse_args()
    selectROI(args.vcf, args.bam, args.bed, args.prefix)

