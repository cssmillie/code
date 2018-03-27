import argparse
import ssub3 as ssub

# get input arguments
parser = argparse.ArgumentParser()
parser.add_argument('-i', help='run directory')
parser.add_argument('-s', help='SampleSheet.csv')
parser.add_argument('-m', help='memory', type=int, default=150)
parser.add_argument('-u', help='usernames', default='')
args = parser.parse_args()

# strip trailing /
rundir = args.i.rstrip('/')

# run bcl2fastq
cmd = 'bcl2fastq --runfolder-dir %s --sample-sheet %s --output-dir ./Data --mask-short-adapter-reads 10 --minimum-trimmed-read-length 10 --no-lane-splitting' %(rundir, args.s)

# submit to cluster
Submitter = ssub.Submitter()
Submitter.H = 'use .bcl2fastq2-2.17.1.14'
Submitter.m = args.m
Submitter.o = 'demult'
Submitter.t = '24:00:00'
Submitter.users = args.u.split(',')
Submitter.submit([cmd])
