#!/usr/bin/python

import sys
import os


def process(lines=None):
    ks = ['name', 'sequence', 'optional', 'quality']
    return {k: v for k, v in zip(ks, lines)}

try:
    fn = sys.argv[1]
except IndexError as ie:
    raise SystemError("Error: Specify file name\n")

if not os.path.exists(fn):
    raise SystemError("Error: File does not exist\n")

record = {}

n = 4
with open(fn, 'r') as fh:
    lines = []
    for line in fh:
        lines.append(line.rstrip())
        if len(lines) == n:
            record = process(lines)
            lines = []

try:
    fq = sys.argv[2]
except IndexError as ie:
    raise SystemError("Error: Specify file name\n")


try:
    t = int(sys.argv[3])
except IndexError as ie:
    raise SystemError("Error: Cutting length not given\n")


subsequences = [record['sequence'][i:i+t] for i in range(0, len(record['sequence']), t)]
subqualities = [record['quality'][i:i+t] for i in range(0, len(record['quality']), t)]

f = open(fq, "w")
increment = 0
for seq,qual in zip(subsequences,subqualities):
    increment = increment +1
    f.write(record['name'] + "-" +str(increment) + "\n")
    f.write(seq+"\n")
    f.write("+\n")
    f.write(qual+"\n")






