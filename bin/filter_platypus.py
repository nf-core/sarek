#!/usr/bin/env python
############################################################
# filter_platypus.py Daniel Nichol 01/05/2018
#
# Filters the vcf file produced by platypus
#
# The filtering criteria are:
#   1) The FILTER tag for the call is in filterList (below)
#   2) The variant is not a known GL variant, or aligned to the decoy genome.
#   3) A genotype is called for all samples (i.e. no sample assigned './.')
#   4) The genotype phred scores (GQ) are >=10 for all samples
#   5) The number of reads covering the variant site is >=10 in all samples.
#   6) The normal sample has no reads exhibing the variant
#   7) At least one tumour samples has >=3 reads for the variant.
#
############################################################
import sys
import os

#####################################################################
# Filter list - the list of FILTER entries that are acceptable.
# Amend this list to include/exclude specific FILTER values.
#####################################################################
filterList =   ['PASS', 'alleleBias', 'Q20', 'Q20;alleleBias','QD', 'Q20;QD', 'QD;alleleBias',
                'Q20;QD;alleleBias', 'SC', 'SC;Q20', 'SC;alleleBias', 'SC;Q20;alleleBias', 'SC;QD',
                'SC;Q20;QD', 'SC;QD;alleleBias', 'SC;Q20;QD;alleleBias', 'HapScore', 'Q20;HapScore',
                'HapScore;alleleBias', 'Q20;HapScore;alleleBias', 'QD;HapScore', 'Q20;HapScore;QD',
                'QD;HapScore;alleleBias', 'Q20;HapScore;QD;alleleBias', 'SC;HapScore', 'SC;Q20;HapScore',
                'SC;HapScore;alleleBias', 'SC;Q20;HapScore;alleleBias', 'SC;HapScore;QD', 'SC;Q20;HapScore;QD',
                'SC;HapScore;QD;alleleBias', 'SC;Q20;HapScore;QD;alleleBias']


#####################################################################
# Parse the arguments
#####################################################################
if len(sys.argv)!=3:
    print("python filter_platypus.py <input> <normal_name>")
    exit()
elif len(sys.argv) == 3:
    platypus_file = sys.argv[1]
    normal_name = sys.argv[2]

#####################################################################
# Filter the platypus file
#####################################################################
with open(platypus_file[0:-4]+"_filtered.vcf",'w') as platypus_pass:
    with open(platypus_file[0:-4]+"_removed.vcf",'w') as platypus_nopass: 
	with open(platypus_file,'r') as platcalls:

            #Copy the header of the platypus file.
            line = ''
            while line[0:6]!="#CHROM":
                line = platcalls.readline()
                platypus_pass.write(line)
                platypus_nopass.write(line)
         
            #Identify the normal sample location in the headers
            header=line[:-1].split('\t')[9:]
            samples= len(header)
            normIx = header.index(normal_name)
    
            line = platcalls.readline()
            while line:
                record = line[:-1].split('\t')
                passed = False

                #Split the (colon (:) separated) individual sample information for each sample.
                for k in range(len(record[:-samples]), len(record)):
                    record[k] = record[k].split(':')
                    # If multiple alts exist, we take NR, NV to be the maximum amongst them
                    if ',' in record[k][4]:
                        record[k][4] = max(record[k][4].split(','))
                        record[k][5] = max(record[k][5].split(','))
    
                allSamples = record[-samples:]
                tumSamples = [allSamples[i] for i in range(len(allSamples)) if i!=normIx]
                normSample = allSamples[normIx]
    
                #Filters: FILTER passed and not decoy/germline
                if (record[6] in filterList) and (record[0]!='hs37d5') and (record[0][0:2]!='GL'):
                    GQsPass = all([float(el[3]) >= 10 for el in allSamples])
                    NRsPass = all([float(el[4]) >= 5 for el in allSamples])
    
                    # GQs (genotype phred scores) all >= 10, and
                    # NRs (number of reads at site) all >= 0
                    if GQsPass and NRsPass:
                        allGTs = [el[0] for el in allSamples]
                        tumGTs = [el[0] for el in tumSamples]
                        normGT = normSample[0]
    
                        # If the genotype for all samples exist:
                        if ("./." not in allGTs) and (normGT=="0/0") and (tumGTs.count("0/0")!=(samples-1)):
    
                            # No variant reads for normal and
                            # >=3 variant reads for at least one tumour sample
                            normNVpass = (int(normSample[5]) == 0)
                            tumNVpass = any([float(el[5]) >= 0 for el in tumSamples])
    
                            if normNVpass and tumNVpass:
                                passed = True
    
                if passed:
                    platypus_pass.write(line)
                else:
                    platypus_nopass.write(line)

                line = platcalls.readline()
