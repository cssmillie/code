import glob, re, sys

regex = sys.argv[1]

for fn in glob.glob(regex):
    for line in open(fn):
        line = line.rstrip()
        if line.startswith('>'):
            line = '>' + fn + '.' + line[1:]
        print(line)
