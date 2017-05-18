#!/usr/bin/env python3
"""
Parses a trace.txt file in which RunHaplotypecaller processes are recorded.
Outputs a list of intervals sorted by runtime, longest duration first.

Run like this:

  ./sort_intervals.py < trace.txt > sorted.list
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


times = []
for line in sys.stdin:
    fields = line.split('\t')
    name = fields[3]
    if not name.startswith('RunHaplotypecaller'):
        continue
    duration = fields[14]
    s = fields[3]
    interval = s[s.find('chr'):-1].replace('_', ':')
    times.append((parse_duration(duration), interval))
times.sort(reverse=True)
for duration, interval in times:
    print(interval)
