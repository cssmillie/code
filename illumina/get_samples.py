import argparse, glob, os, re

# This automatically constructs a samples.txt file for the RNA-Seq pipeline
# It searches for demultiplexed forward and reverse FASTQ files
# Then uses the groups from the "--regex" argument to get the sample names

# Usage example:
# python ~/code/illumina/get_samples.py --data "170602_NB501164_0451_AHYLFMBGXY/Data" --regex '.*GFP_(.*)_sNUCSeq.*(_S\d+).*'


# Get input arguments
parser = argparse.ArgumentParser()
parser.add_argument('--data', help='Data folder (demultiplex output)')
parser.add_argument('--regex', help='Regex to build sample names')
parser.add_argument('--join', help='String to join groups with', default='')
args = parser.parse_args()

# Get forward and reverse reads
fwd = sorted(glob.glob('%s/*/*/*R1*fastq.gz' %(args.data)))
rev = sorted(glob.glob('%s/*/*/*R2*fastq.gz' %(args.data)))

# Get sample names
for i in range(len(fwd)):
    fi = os.path.abspath(fwd[i])
    ri = os.path.abspath(rev[i])
    assert re.sub('R1', '', fi) == re.sub('R2', '', ri)
    name = args.join.join(re.search(args.regex, fi).groups())
    print '%s\t%s\t%s' %(name, fi, ri)
