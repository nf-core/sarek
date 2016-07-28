from Bio.SeqIO.QualityIO import FastqGeneralIterator
import sys 
ids = set(x[:-1] for x in open(sys.argv[1]))
for name, seq, qual in FastqGeneralIterator(sys.stdin):
    if (name in ids): print("@%s\n%s\n+\n%s" % (name, seq, qual))
