#!/usr/bin/env python
#
# The MIT License (MIT)
#
# Copyright (c) 2018 Szilveszter Juhos
# the reduce() code shamelessly copied from
# https://github.com/brentp/interlap
# Copyright (c) 2014 Brent S. Pedersen
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.

import os
import sys
import copy
import argparse
import subprocess

def parse_args():
    parser = argparse.ArgumentParser(description='Selecting Region Of Interest from BAM files')
    parser.add_argument('-v', help='List of VCF files to process like -v tumour.vcf,normal.vcf', dest="vcfs")
    parser.add_argument('-a', help='The BAM files (indexed) to read from like -b tumour.bam,normal.bam', required=True, dest="bams")
    parser.add_argument('-b', help='The BED file with gene coordinates', dest="beds")
    parser.add_argument('-p', help='The prefix of the target BAM file to write to (postfix is the input filename)', dest="prefix", default="ROI.")
    parser.add_argument('-t', help='Number of threads used by samtools [8]', dest="threads", default=8, type=int)
    parser.add_argument('-w', help='Region width for VCF loci [200]', dest="width", default=200, type=int)
    parser.add_argument('-r', help='Readlength for region padding [151]', dest="pad", default=151, type=int)
    parser.add_argument('-c', help='Number of variants in a chunk [150]', dest="maxRegions", default=150, type=int)

    return parser.parse_args()

def reduce(args):
    """
    >>> reduce([(2, 4), (4, 9)])
    [(2, 4), (4, 9)]

    >>> reduce([(2, 6), (4, 10)])
    [(2, 10)]
    """
    if len(args) < 2: return args
    args.sort()
    ret = [args[0]]
    for next_i, (s, e) in enumerate(args, start=1):
        if next_i == len(args):
            ret[-1] = ret[-1][0], max(ret[-1][1], e)
            break

        ns, ne = args[next_i]
        if e > ns or ret[-1][1] > ns:
            ret[-1] = ret[-1][0], max(e, ne, ret[-1][1])
        else:
            ret.append((ns, ne))
    return ret

class ROISelector:
    """Class to selet region of interest (ROI) from a set of BAM files. ROI is defined by input VCFs and BEDs""" 

    def __init__(self, threads, maxRegions, prefix, width, pad):
        # this is the set where loci are stored
        # for every chromosome there is a list in the dict
        # VCF records are int-s in the list,
        # BED entries are (start,end) tuples
        self.callDict = dict()
        # number of samtools threads
        self.threads = threads
        # number of variants in a chunk
        self.maxRegions = maxRegions
        # prefix of the ROI BAM file
        self.prefix = prefix
        # width of region for VCF records
        self.pad = pad
        # add readlength padding
        self.width = width + pad

    def selectROI(self,vcfs,bams,beds):
        if not vcfs and not beds:
            print("At least a single VCF (-v) or BED (-b) file is expected")
            sys.exit()

        # first we MUST process BED intervals, since we want to ignore all the 
        # variants that are already in any of the intervals
        bedRecordCount = 0
        if beds:                 # if there are any BEDs 
            print("Processing BED files: ")
            bedFiles = beds.split(",")
            for bedLoci in bedFiles:
                print("\t" + bedLoci)
                self.addBEDRecords(bedLoci)
            for chrom in self.callDict:
                bedRecordCount = bedRecordCount + len(self.callDict[chrom])
            print("Collected " + str(bedRecordCount) + " entries from BEDs")

        if vcfs:                # if there are any VCFs given
            print("Processing VCF files: ")
            vcfFiles = vcfs.split(",")
            for vcfCalls in vcfFiles:
                print("\t" + vcfCalls)
                self.addVCFRecords(vcfCalls)
            vcfRecordCount = 0
            # report the number of records added
            for chrom in self.callDict:
                vcfRecordCount = vcfRecordCount + len(self.callDict[chrom])
            print("Collected " + str(vcfRecordCount-bedRecordCount) + " entries from VCFs")
        # now collapse overlapping intervals
        for chrom in self.callDict:
            self.callDict[chrom] = reduce( self.callDict[chrom] )
            # we still have to collapse intervals that are too close to each other 
            # if they are closer than the readlength, there can be cases when a read is picked 
            # in both interval, so we have to merge them
            self.mergeCloseIntervals(chrom)

        # save intervals to a BAM
        print("Saving reads from BAMs:")
        self.saveCallsToBAM(bams)

    def mergeCloseIntervals(self,chrom):
        # we are assuming the intervals are ordered (they should to be)
        mergedChrom = []
        start = (0,0)
        for interv in self.callDict[chrom]:
            # (start[0], start[1])
            #                          (interv[0], interv[1])
            # |------------------|     |--------------------|     
            if interv[0] - start[1] <= self.pad*2:
                # merge the two intervals
                start = (start[0],interv[1])
            elif start[1]!= 0:  # if its not the very first one and have enough space between the intervals
                # add to the new interval set
                mergedChrom.append(start)
                start = interv
            else: # most of the intervals should fall here
                start = interv
        self.callDict[chrom] = mergedChrom

    def writeChunk(self,bam,chrom,regions):
        # from the regions list make a list of lists - sublists have 50 entries max
        sublists = [regions[i:i+self.maxRegions] for i  in range(0, len(regions), self.maxRegions)]
        regionCount = 0
        filesToMerge = ""
        for rlist in sublists:
            chunkBAM = chrom +"_"+ str(regionCount) + "_tmpROI.bam "
            regionCount = regionCount + 1
            filesToMerge = filesToMerge + chunkBAM
            loci = ""
            for locus in rlist:
                loci = loci + locus + " "
            subprocess.call("samtools view -@" + str(self.threads) + " " + bam + " -bo " + chunkBAM + loci, shell=True)
            print("Wrote chunk for " + chrom + " to " + chunkBAM)
        return filesToMerge

    def saveCallsToBAM(self,bams):
        bamFiles = bams.split(",")
        for bam in bamFiles:
            print("\t" + bam)
            # getting reads one-by-one takes aages. We are going to try to use htslib/samtools
            # by making a command line to get ROI for each chromosome, and 
            # merge/sort the BAMs later
            filesToMerge = ""
            for chrom in self.callDict:
                regionCount = 0;
                regions = []
                for gcoord in self.callDict[chrom]:
                    if type(gcoord) is tuple:       # BED record in tuple
                        start = gcoord[0]
                        end = gcoord[1]
                    else:                   # VCF record
                        start = gcoord - self.width
                        end = gcoord + self.width
                    regions.append(chrom + ":" + str(start) + "-" + str(end))
                # let's write the file for the chromosome
                filesToMerge = filesToMerge + self.writeChunk(bam,chrom,regions)
            print("Merging ...")
            subprocess.call("samtools merge -@"+str(self.threads)+" -f merged.roi.bam " + filesToMerge, shell=True)
            target = os.path.dirname(bam) 
            if len(target) > 0:
                target = target + "/"
            target = target + self.prefix + os.path.basename(bam)
            print("Sorting ...")
            subprocess.call("samtools sort -@" + str(self.threads) + " -o " + target + " merged.roi.bam", shell=True)
            print("Sorted file written to " + target)
            print("Indexing ...")
            subprocess.call("samtools index -@" + str(self.threads) + " " + target, shell=True)

    def addVCFRecords(self,VCFFile):
        """
        Store VCF record CHROM,POS in a dict 
        """
        VCFlines = [line.rstrip() for line in open(VCFFile,'rt')]
        for variant in VCFlines:
            if not variant.startswith("#"):      # ignore header
                vcfCols = variant.split()
                chrom = vcfCols[0]
                locus = int(vcfCols[1]) 
                if chrom not in self.callDict:   # if chromosome not in dictionary yet
                    self.callDict[chrom] = []
                self.callDict[chrom].append( (locus-self.width, locus+self.width) )

    def addBEDRecords(self,BEDFile):
        """
        Store BED record CHROM,START,END in a set
        """
        BEDlines = [line.rstrip() for line in open(BEDFile,'rt')]
        for locus in BEDlines:
            BEDrecord = locus.split()
            chrom = BEDrecord[0] 
            if chrom not in self.callDict:
                self.callDict[chrom] = []
            self.callDict[chrom].append( (int(BEDrecord[1]),int(BEDrecord[2])) )

######################## end of class ROISelector ######################################


if __name__ == "__main__":
    args = parse_args()
    rois = ROISelector(args.threads, args.maxRegions, args.prefix, args.width, args.pad)
    rois.selectROI(args.vcfs,args.bams,args.beds)

