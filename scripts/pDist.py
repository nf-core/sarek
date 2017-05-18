#!/usr/bin/env python
import argparse

def parse_args():
    parser = argparse.ArgumentParser(description='Checking call density for a sorted VCF')
    parser.add_argument('-v', help='A sorted VCF file to check',required=True, dest="vcf")
    parser.add_argument('-s', help='Step size (in million basepairs [5000000])', required=False, dest="step", default=5000000)
    parser.add_argument('--nochrom', help='Are we expecting chr at the beginning of the chromosome string? [True]', action="store_true")
    return parser.parse_args()


def printVCFdensity(vcf,step,nochrom):
    with open(vcf,"r") as fh:
        currentStep = step      # current coordinates max value
        varBin = 0              # number or variants in the current bin
        currentChromosome = ""  # ID of the current chromosome
        # put chr into chromosome names if expected
        chromosomes = map((lambda x: str(x)),range(1,23) + ['X','Y']) if nochrom else map((lambda x: 'chr'+str(x)),range(1,23) + ['X','Y'])

        for line in fh:
            if not line.startswith("#"):
                parts = line.split()
                chrom = parts[0]
                coord = int(parts[1])
                if chrom in chromosomes:
                    if chrom != currentChromosome:
                        if varBin < 100 and currentChromosome != "":
                            print currentStep, varBin
                        currentChromosome = chrom
                        #printBar(currentStep, varBin)
                        # restart counting
                        currentStep = step
                        varBin = 0
                        print "processing " + currentChromosome
                    if currentStep < coord:
                        if varBin < 100:
                            print currentStep, varBin
                        #printBar(currentStep, varBin)
                        currentStep += step
                        varBin = 1
                    else:   # the actual coordinate is smaller or equal to the current step max
                        varBin += 1

def printBar(s,vb):
    bar = "#"
    for hm in range(1,int((float(vb)/100))):
        bar += "#"
    print s,bar

if __name__ == "__main__":
    args = parse_args()
    printVCFdensity(args.vcf,args.step,args.nochrom)
