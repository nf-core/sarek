#!/usr/bin/env python
import click

# consider chromosome names starting with chr
chromosomes = map((lambda x: 'chr'+str(x)),range(1,23) + ['X','Y'])
# we are going to have a dictionary,  where the first keys are chromosome names, the second are 
# the ranges (max values as we are expecting sorted VCFs)
# { 'chr1': {10000000: 123, 20000000: 345, ...}, 'chr2': {}


@click.command(context_settings = dict( help_option_names = ['-h', '--help'] ))
@click.option('--vcf','-v', type=str, help='A sorted VCF file to check')
@click.option('--step','-s', type=int, default=5000000, help='Step size (in million basepairs [5000000])')
def printVCFdensity(vcf,step):
    with open(vcf,"r") as fh:
        currentStep = step      # current coordinates max value
        varBin = 0              # number or variants in the current bin
        currentChromosome = ""  # ID of the current chromosome
        for line in fh:
            if not line.startswith("#"):
                parts = line.split()
                chrom = parts[0]
                coord = int(parts[1])
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
    printVCFdensity()
