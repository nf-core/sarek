#!/usr/bin/env python
import os
import sys
import copy
import argparse
import subprocess

def parse_args():
    parser = argparse.ArgumentParser(description='Selecting Region Of Interest from BAM files')
    parser.add_argument('-v', help='List of VCF files to process like -v tumour.vcf,normal.vcf',required=False, dest="vcfs")
    parser.add_argument('-a', help='The BAM files (indexed) to read from like -b tumour.bam,normal.bam', required=True, dest="bams")
    parser.add_argument('-b', help='The BED file with gene coordinates', required=False, dest="beds")
    parser.add_argument('-p', help='The prefix of the target BAM file to write to (postfix is the input filename)', required=False, dest="prefix", default="ROI.")
    parser.add_argument('-t', help='Number of threads used by samtools [8]', required=False, dest="threads", default=8)
    parser.add_argument('-w', help='Region width for VCF loci [200]', required=False, dest="width", default=200)
    parser.add_argument('-c', help='Number of variants in a chunk [150]', required=False, dest="maxRegions", default=150)

    return parser.parse_args()


class bedChunk:
    """Small class to define a genomic locus with chromosome and start and end"""
    def __init__(self):
        self.loci = dict()

class ROISelector:
    """Class to selet region of interest (ROI) from a set of BAM files. ROI is defined by input VCFs and BEDs""" 

    def __init__(self, threads, maxRegions, prefix, width):
        # this is the set where loci are stored
        # for every chromosome there is a list in the dict
        # VCF records are int-s in the list,
        # BED entries are (start,end) tuples
        self.callDict = dict()
        # number of samtools threads
        self.threads = int(threads)
        # number of variants in a chunk
        self.maxRegions = int(maxRegions)
        # prefix of the ROI BAM file
        self.prefix = prefix
        # width of region for VCF records
        self.width = int(width)
        # we will make a deep copy of the callDict
        # to save the BED intervals into a separate entry
        self.callDictCopy = None

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
        # we make this copy to process intervals faster
        self.callDictCopy = copy.deepcopy(self.callDict)
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
        # save intervals to a BAM
        print("Saving reads from BAMs:")
        self.saveCallsToBAM(bams)

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
            target = os.path.dirname(bam) + "/" + self.prefix + os.path.basename(bam)
            print("Sorting ...")
            subprocess.call("samtools sort -@" + str(self.threads) + " -o " + target + " merged.roi.bam", shell=True)
            print("Sorted file written to " + target)
            print("Indexing ...")
            subprocess.call("samtools index -@" + str(self.threads) + " " + target, shell=True)

    def isInBEDIntervals(self,chrom,locus):
        """"We are checking whether the locus is already covered by a BED interval"""
        # the copy contains only the BED intervals
        # and we are hoping there are only few intervals, since for many intervals this search is not fast enough
        BEDloci = self.callDictCopy[chrom]
        for region in BEDloci:      # these are tuples
            if  region[0] <= locus and locus <= region[1]:
                #print("Locus "+chrom +":"+str(locus) + " is in region " + chrom + ":"+ str(region[0]) + "-" + str(region[1]))
                return True
        return False

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
                if not self.isInBEDIntervals(chrom,locus):
                    self.callDict[chrom].append(locus)

    def addBEDRecords(self,BEDFile):
        """
        Store BED record CHROM,START,END in a set
        """
        BEDlines = [line.rstrip() for line in open(BEDFile,'rt')]
        for locus in BEDlines:
            BEDrecord = locus.split()
            if BEDrecord[0] not in self.callDict:
                self.callDict[BEDrecord[0]] = []
            self.callDict[BEDrecord[0]].append( (int(BEDrecord[1]),int(BEDrecord[2])) )

######################## end of class ROISelector ######################################


if __name__ == "__main__":
    args = parse_args()
    rois = ROISelector(args.threads, args.maxRegions, args.prefix, args.width)
    rois.selectROI(args.vcfs,args.bams,args.beds)

