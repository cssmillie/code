import argparse, ssub

# get input arguments
parser = argparse.ArgumentParser()
parser.add_argument('--dir', help='run directory')
args = parser.parse_args()

# strip trailing /
rundir = args.dir.rstrip('/')

# run bcl2fastq
cmd = 'bcl2fastq --runfolder-dir %s --output-dir %s/Data --mask-short-adapter-reads 10 --minimum-trimmed-read-length 10 --no-lane-splitting' %(rundir, rundir)

# submit to cluster
Submitter = ssub.Submitter()
Submitter.H = 'use .bcl2fastq2-2.17.1.14'
Submitter.m = 150
Submitter.o = 'demult'
Submitter.submit(cmd)
