import argparse, glob, numpy, os, re

# input arguments
parser = argparse.ArgumentParser()
parser.add_argument('--dir', help='demultiplexd sequences folder')
parser.add_argument('--csv', help='samples.csv (lane, sample, index)')
parser.add_argument('--out', help='output folder')
parser.add_argument('--lns', help='symbolic link', default=False, action='store_true')
parser.add_argument('--keep', help='move files in csv', default=False, action='store_true')
args = parser.parse_args()

# map: sample index -> oligonucleotides
# e.g. SI-3A-A1 -> AAACGGCG,CCTACCAT,GGCGTTTC,TTGTAAGA
m = {}
for line in open('/home/unix/csmillie/code/10X/sample_indices.txt'):
    line = line.rstrip().split(',')
    m[line[0]] = line[1:]

# split fastqs into groups by lane
def split_by_lane(fastqs, lane):
    regex = '_lane-0*%s' %(re.sub('NA', '[0-9]*', lane))
    keep = []
    remove = []
    for fastq in fastqs:
        if re.search(regex, fastq):
            keep.append(fastq)
        else:
            remove.append(fastq)
    u = sum([os.path.getsize(fn) for fn in keep])
    v = sum([os.path.getsize(fn) for fn in remove])
    print 'keep: %d, remove: %d' %(len(keep), len(remove))
    print 'lane %s uses %s%% of data' %(lane, 100.*u/(u+v))
    return [keep, remove]

# search for files with lane and index
def get_fastqs(lane, index):
    fastqs = []
    for oligo in m[index]:
        i1, I1 = split_by_lane(glob.glob('%s/outs/fastq_path/read-I1_si-%s_lane*fastq.gz' %(args.dir, oligo)), lane)
        i2, I2 = split_by_lane(glob.glob('%s/outs/fastq_path/read-I2_si-%s_lane*fastq.gz' %(args.dir, oligo)), lane)
        ra, RA = split_by_lane(glob.glob('%s/outs/fastq_path/read-RA_si-%s_lane*fastq.gz' %(args.dir, oligo)), lane)
        assert len(i1) == len(i2) == len(ra)
        assert len(I1) == len(I2) == len(RA)
        fastqs += (i1 + i2 + ra)
    return fastqs

# main routine
keep = []
for line in open(args.csv):
    if line.startswith('Lane'):
        continue
    lane, sample, index = line.rstrip().split(',')
    fastqs = get_fastqs(lane, index)
    keep += fastqs

remove = [fn for fn in glob.glob('%s/outs/fastq_path/*.fastq.gz' %(args.dir)) if fn not in keep]

print keep

print 'Keeping %d files (total size = %d)' %(len(keep), sum([os.path.getsize(fn) for fn in keep]))
print 'Removing %d files (total size = %d)' %(len(remove), sum([os.path.getsize(fn) for fn in remove]))

# moving files! be extra careful!
if not os.path.exists(args.out):
    os.mkdir(args.out)

# get files to move
if args.keep:
    move = keep
else:
    move = remove

mv = open('%s/mv.sh' %(args.out), 'w')
for fn in move:
    new_fn = '%s' %(os.path.basename(fn))
    if os.path.exists(new_fn):
        quit('Error! destination file %s already exists! Do not clobber data!' %(new_fn))
    if args.lns:
        mv.write('ln -s ../%s ./%s\n' %(fn, new_fn))
    else:
        mv.write('mv ../%s ./%s\n' %(fn, new_fn))
mv.close()
    
