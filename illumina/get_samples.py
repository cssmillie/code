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
parser.add_argument('--cells', help='Number of expected cells', default=None, type=int)
parser.add_argument('--cloud', help='Copy data to google cloud', default=False, action='store_true')
parser.add_argument('--out', help='Output file (samples.txt)', default='')
args = parser.parse_args()

def get_cpath(fn):
    # get path on google cloud
    return '%s%s' %(os.environ['CLOUD'], os.path.abspath(fn))

# Check input arguments
if args.cloud:
    if args.out == '':
        quit('Must specify output file (--out)')

# Get forward and reverse reads
fwd = sorted(glob.glob('%s/*/*/*R1*fastq.gz' %(args.data)))
rev = sorted(glob.glob('%s/*/*/*R2*fastq.gz' %(args.data)))
if args.out:
    out = open(args.out, 'w')
    out.write('Cell,Plate,Read1,Read2\n')
else:
    print('Cell,Plate,Read1,Read2\n')

# Get sample names
for i in range(len(fwd)):
    fi = Fi = os.path.abspath(fwd[i])
    ri = Ri = os.path.abspath(rev[i])
    assert re.sub('R1', '', fi) == re.sub('R2', '', ri)
    if args.cloud:
        Fi = get_cpath(fi)
        Ri = get_cpath(ri)
    name = args.join.join(re.search(args.regex, fi).groups())
    line = '%s,%s,%s,%s' %(name, 'plate', Fi, Ri)
    if args.out:
        out.write(line + '\n')
    else:
        print(line)

if args.cloud:
    out.close()
    
    cmd = 'gsutil -m cp -r %s %s' %('Data', re.sub('Data', '', get_cpath('Data')))
    print(cmd)
    os.system(cmd)
    
    cmd = 'gsutil -m cp %s %s' %(args.out, get_cpath(args.out))
    print(cmd)
    os.system(cmd)
    
    print(get_cpath(args.out))
    print(get_cpath('.'))

