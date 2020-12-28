#!/usr/bin/env python
import sys

# filter blast report

x = {}

for line in sys.stdin:
    line = line.rstrip().split('\t')

    sid = line[0]
    pct = float(line[2])
    evalue = float(line[10])

    x[sid] = x.get(sid, 0) + 1

    if (x[sid] <= 10) and ((pct >= 95) and (evalue <= .0001)):
        print('\t'.join(line))
