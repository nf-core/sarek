#!/usr/bin/env python
######################################################################
# 23/08/2018 - Dan Nichol & Timon Heide
#
# Filters a VCF files to discard any calls that:
#   1. Do not satisfy the filter criteria of Mutect2 (i.e. not PASS)
#   2. Align to a non canonical chromosome (i.e (chr)1-22)
#   3. Do not have total allelic depth (AD) >=10 in both tumour and normal
#   4. Do not satisfy (ALT_F1R2 + ALT_F2R1)>=3 in the tumour sample
######################################################################

import sys
import re

if len(sys.argv) != 3:
    print("Usage: ", sys.argv[1], "<input> <output>")

vcfFileInput = sys.argv[1]
vcfFileOutput = sys.argv[2]


with open(vcfFileInput,'r') as in_vcf:
    with open(vcfFileOutput,'w') as out_vcf:

        # parse the header:
        for header_line in in_vcf:
            out_vcf.write(header_line)

            if header_line.startswith("##normal_sample="):
                normal_id = header_line[len("##normal_sample="):][:-1]

            if header_line.startswith("##tumor_sample="):
                tumour_id = header_line[len("##tumor_sample="):][:-1]

            if header_line.startswith("#CHROM"): # this marks end of header
                row_names = header_line[1:-1].split('\t')
                break

        # parse mutations
        for mutation_line in in_vcf:
            fields = dict(zip(row_names, mutation_line[:-1].split('\t')))
            info = fields["FORMAT"].split(':')
            n_info = dict(zip(info, fields[normal_id].split(':')))
            t_info = dict(zip(info, fields[tumour_id].split(':')))
            n_info["AD"] = [int(e) for e in n_info["AD"].split(',')]
            t_info["AD"] = [int(e) for e in t_info["AD"].split(',')]

            pass_filter = fields["FILTER"] == "." or fields["FILTER"] == "PASS"
            pass_chromo =  bool(re.match("^(chr)?[0-9XY]*$", fields["CHROM"]))

            if pass_chromo \
               and pass_filter \
               and sum(t_info["AD"]) > 10 \
               and sum(n_info["AD"]) > 10 \
               and t_info["AD"][1] >= 3 \
               and n_info["AD"][1] == 0 \
               and (t_info["GT"] != "0/0" and t_info["GT"] != "0|0") \
               and (n_info["GT"] == "0/0" or n_info["GT"] == "0|0"):
                    out_vcf.write(mutation_line)
