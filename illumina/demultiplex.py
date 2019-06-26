import argparse
import ssub3 as ssub

# get input arguments
parser = argparse.ArgumentParser()
parser.add_argument('-i', help='run directory')
parser.add_argument('-s', help='SampleSheet.csv')
parser.add_argument('-m', help='memory', type=int, default=150)
parser.add_argument('-u', help='usernames', default='')
parser.add_argument('-l', help='lanes', default='') # range or comma-separated list
parser.add_argument('-p', help='print', default=False, action='store_true')
args = parser.parse_args()

# strip trailing /
rundir = args.i.rstrip('/')

# run bcl2fastq
cmd = 'bcl2fastq --runfolder-dir %s --sample-sheet %s --output-dir ./Data --mask-short-adapter-reads 10 --minimum-trimmed-read-length 10' %(rundir, args.s)

if args.l:
    cmd = cmd + ' --tiles s_[%s]' %(''.join(args.l))
cmd = cmd + ' --no-lane-splitting'

print(cmd)
if args.p:
    quit()

# submit to cluster
Submitter = ssub.Submitter()
#Submitter.H = 'use .bcl2fastq2-2.19.1.403' # redhat6
#Submitter.H = 'use .bcl2fastq2-2.20.0.422' # redhat6 and redhat7
Submitter.H = 'use .bcl2fastq2-2.20.0.422; use .bcl2fastq2-v2.20.0' # redhat 6 and redhat7
Submitter.m = args.m
Submitter.o = 'demult'
#Submitter.O = 'RedHat6'
Submitter.t = '24:00:00'
Submitter.users = args.u.split(',')
Submitter.submit([cmd])
