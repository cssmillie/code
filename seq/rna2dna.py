import re, sys, util

fn = sys.argv[1]

for record in util.iter_seq(fn):
    record[1] = re.sub('U', 'T', record[1])
    print('\n'.join(record))
