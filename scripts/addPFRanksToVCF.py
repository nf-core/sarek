#!/usr/bin/python3.6
"""
Pathfindr (https://github.com/NBISweden/pathfindr) stores ranked annotation in a CSV file. 
This script extracts these ranks scores and chromosome coordinates, and adds the scores
as annotations to the input VCF.

Inputs: 
  -v  --vcf     the VCF to annotate
  -c  --csv     the CSV with rank scores
  -f --family   family ID (normal sample name)

Results are written to STDOUT
"""

import click
import re


class addPFRanksToVCF:

  def __init__(self,familyID):
    self.ranks = {}             # dict with a (chr,coord) as key, and rank as value
    self.rank_idx = 0           # column for rank_score (first line have to be parsed to find it)
    self.chr_idx = 0            # ditto, but for chromosome
    self.start_idx = 0          # and for the start coordinate
    self.family_id = familyID   # sample name - usually the normal

    self.header_line = re.compile("^#")
    self.chrom_line = re.compile("^#CHROM")

  def readRanks(self,csv):
    with open(csv,'r') as infile:
      for line in infile:
        # we have to parse the very first line to find the columns for values
        if 'rank_score' in line:
          self.findColumns(line)
        else:
          cols = line.split(",")
          score = cols[self.rank_idx]
          chrom = cols[self.chr_idx]
          coord = cols[self.start_idx]
          self.ranks[(chrom,coord)] = score

  def printVCF(self,vcf):
    with open(vcf,'r') as vcffile:
      self.printHeader(vcffile)
      for line in vcffile:
        # TODO skip in printHeader()
        if not self.header_line.match(line):
          vcfcols = line.split()
          rank_score = 0;
          if (vcfcols[0],vcfcols[1]) in self.ranks:
            rank_score = self.ranks[(vcfcols[0],vcfcols[1])]
          self.addRankToLine(vcfcols,rank_score)
 
  def printHeader(self,vcffile):
    for line in vcffile:
      if self.chrom_line.match(line):
        print("##INFO=<ID=RankScore,Number=.,Type=String,Description=\"The (somatic) rank score for this variant calculated by https://github.com/NBISweden/pathfindr. family_id:rank_score.\">")
        line = line.replace("NORMAL",str(self.family_id))
        print(line,end="")
        return
      if self.header_line.match(line):
        print(line,end="")
    
  def addRankToLine(self,cols,rank):
    """
    We assume VCF lines are missing the RankScore annotation, so we are giving one.
    It has to be something like ;RankScore=FAMILY:2 added to the INFO (8th) column
    """
    line = ""
    ann = ";RankScore=" + str(self.family_id) + ":" + str(rank)
    for c in cols:
      line = line + c
      if cols.index(c) == 7:  # 0-based index
        line = line + ann
      if cols.index(c) < len(cols):
        line = line + "\t"
    print(line)

  def findColumns(self,line):
    cols = line.split(",")
    self.rank_idx = cols.index("rank_score")
    self.chr_idx = cols.index("chr")
    self.start_idx = cols.index("start")


# This is the surrogate for main(): everything happens here

@click.command(context_settings = dict( help_option_names = ['-h', '--help'] ))
@click.option('--vcf',      '-v', type=str, help='VCF file to annotate', required=True)
@click.option('--csv',      '-c', type=str, help='CSV file with ranked scores', required=True)
@click.option('--family',   '-f', type=str, help='family ID as expected in scout yaml file (usually sample name of the normal)', required=True)

def annotateVCF(vcf,csv,family):
  ranker = addPFRanksToVCF(family)

  # make a dict with coords and ranks
  ranker.readRanks(csv)
  ranker.printVCF(vcf)

if __name__ == "__main__":
  annotateVCF()
