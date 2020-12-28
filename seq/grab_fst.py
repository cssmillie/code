#!/usr/bin/python

import re, sys

regex = sys.argv[1]

flag = 0
for line in sys.stdin:
    line = line.rstrip()
    if line.startswith('>'):
        if re.search(regex, line):
            flag = 1
        else:
            flag = 0
    if flag == 1:
        print(line)
