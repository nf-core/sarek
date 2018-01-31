#!/usr/bin/env python
import sys
import os
import click
import pysam

@click.command(context_settings = dict( help_option_names = ['-h', '--help'] ))
@click.option('--vcf',  '-v', type=str, help='List of VCF files to process like -v tumour.vcf,normal.vcf', required=False)
@click.option('--bam',  '-b', type=str, help='The BAM files (indexed) to read from like -b tumour.bam,normal.bam', required=True)
@click.option('--bed',  '-e', type=str, help='The BED file with gene coordinates', required=False)
@click.option('--out',  '-o', type=str, help='The target BAM file to write to', required=True)
def selectROI(vcf,bam,bed,out):
    vcfFiles = vcf.split(",")
    bamFiles = bam.split(",")
    callSet = set()
    print("Processing VCF files: ")
    for vcfCalls in vcfFiles:
        print("\t" + vcfCalls)
        addRecords(vcfCalls,callSet)
    print("Number of collected records " + str(len(callSet)))
    print("Saving reads from BAMs:")
    saveVCFCallsToBAM(callSet,bamFiles)

def saveVCFCallsToBAM(callSet,bamFiles):
    for bam in bamFiles:
        print("\t" + bam)
        bamfile = pysam.AlignmentFile(bam, "rb")                # the input file
        readSet = set()                                         # this set is to store read IDs: when processing overlapping regions, reads will be written only once
        outBAM = pysam.AlignmentFile("tmpROI.bam","wb",template=bamfile)  # the output file
        for call in sorted(callSet):
            chromosome = call[0]
            start = call[1]-500
            end = call[1]+500
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

def addRecords(callsInFile,callSet):
    """
    Store VCF record CROM,POS in a set
    """

    vcf_in = pysam.VariantFile(callsInFile)
    for record in vcf_in.fetch():
        callSet.add( (record.chrom, record.pos) )

if __name__ == "__main__":
    selectROI()

