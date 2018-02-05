#!/usr/bin/env python
import os
import argparse
import pysam

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
    callSet = set()

    print("Processing VCF files: ")
    vcfFiles = vcf.split(",")
    for vcfCalls in vcfFiles:
        print("\t" + vcfCalls)
        addVCFRecords(vcfCalls,callSet)
    print("Collected " + str(len(callSet)) + " records from VCFs")

    print("Processing BED files: ")
    bedFiles = bed.split(",")
    for bedLoci in bedFiles:
        print("\t" + bedLoci)
        addBEDRecords(bedLoci,callSet)
    print("Collected " + str(len(callSet)) + " records from BEDs")


    bamFiles = bam.split(",")
    print("Saving reads from BAMs:")
    saveCallsToBAM(callSet,bamFiles)

def saveCallsToBAM(callSet,bamFiles):
    for bam in bamFiles:
        print("\t" + bam)
        bamfile = pysam.AlignmentFile(bam, "rb")                # the input file
        readSet = set()                                         # this set is to store read IDs: when processing overlapping regions, reads will be written only once
        outBAM = pysam.AlignmentFile("tmpROI.bam","wb",template=bamfile)  # the output file
        for call in sorted(callSet):
            chromosome = call[0]
            if len(call) == 2:
                # this is a VCF call
                start = call[1]-500
                end = call[1]+500
            else:
                # BED locus
                start = call[1]
                end = call[2]
            # store reads:
            for read in bamfile.fetch(chromosome, start, end):
                read_id = (read.reference_name, read.reference_start, read.query_name)
                if read_id not in readSet:                      # write out read only if not written already
                    readSet.add(read_id)
                    outBAM.write(read)
        outBAM.close()

        # once the filtered reads are written, we have to sort the files
        bn = os.path.basename(bam)                              # we are creating the output file here
        dn = os.path.dirname(bam)
        fn = dn + "/ROI."+bn
        print("Writing intervals to " + fn)
        pysam.sort("-o", fn, "tmpROI.bam") 

def addVCFRecords(VCFFile,callSet):
    """
    Store VCF record CHROM,POS in a set
    """

    vcf_in = pysam.VariantFile(VCFFile)
    for record in vcf_in.fetch():
        callSet.add( (record.chrom, record.pos) )

def addBEDRecords(BEDFile,callSet):
    """
    Store BED record CHROM,START,END in a set
    """

    BEDlines = [line.rstrip() for line in open(BEDFile,'rt')]
    for locus in BEDlines:
        BEDrecord = locus.split()
        callSet.add( (BEDrecord[0],int(BEDrecord[1]),int(BEDrecord[2])) )


if __name__ == "__main__":
    args = parse_args()
    selectROI(args.vcf, args.bam, args.bed, args.prefix)

