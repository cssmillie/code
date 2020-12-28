import sys, util

if len(sys.argv) > 1:
    fn = sys.argv[1]
else:
    fn = sys.stdin

for record in util.iter_fsq(fn):
    print('>%s\n%s' %(record[0][1:], record[1]))
