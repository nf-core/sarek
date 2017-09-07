#!/usr/bin/env python3
"""
Annotates a BED file with RunHaplotypecaller run times

Runtimes come from a trace.txt file.

Run like this:

  python3 runtimeannotate.py intervals.bed trace.txt > annotated.bed
"""
import re
import sys
from itertools import islice


duration_regex = re.compile('((?P<hours>[0-9]*)h($| ))?((?P<minutes>[0-9]*)m($| ))?((?P<seconds>[0-9]*(.[0-9]*)?)s)?')


def parse_duration(s):
    m = duration_regex.match(s)
    if not m:
        return None
    t = 0.
    if m.group('hours'):
        t += 3600 * int(m.group('hours'))
    if m.group('minutes'):
        t += 60 * int(m.group('minutes'))
    if m.group('seconds'):
        t += float(m.group('seconds'))
    return t


if len(sys.argv) != 3:
	print('Pass BED and trace.txt as parameters', file=sys.stderr)
	sys.exit(1)

runtimes = dict()  # Maps (chrom, start, end) to runtime
with open(sys.argv[2]) as trace:
	for line in trace:
		fields = line.split('\t')
		name = fields[3]
		if not name.startswith('RunHaplotypecaller'):
			continue
		duration = fields[14]
		s = fields[3]
		chrom, coords = s[s.find('chr'):-1].split('_')
		start, end = coords.split('-')
		start = int(start) - 1
		end = int(end)
		runtimes[(chrom, start, end)] = parse_duration(duration)

with open(sys.argv[1]) as bed:
	for line in bed:
		fields = line.split('\t')
		chrom = fields[0]
		start = int(fields[1])
		end = int(fields[2])
		d = runtimes.get((chrom, start, end))
		if d is not None:
			#print('{:10.1f} {:9d}'.format((end - start) / d, end-start), chrom, start, end, 'seconds', d, sep='\t')
			print(chrom, start, end, 'seconds', d, sep='\t')
		else:
			print(chrom, start, end)
