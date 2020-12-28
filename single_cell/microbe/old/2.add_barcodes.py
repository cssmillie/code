import argparse, os, re, util

# Get input arguments
parser = argparse.ArgumentParser()
parser.add_argument('--sid', help='Sample ID')
parser.add_argument('--out', help='Output file')
args = parser.parse_args()

# Read Sample Info file
sites = {'EpiA':'Epi_A', 'EpiB':'Epi_B', 'LPA':'LP_A', 'LPB':'LP_B'}
info = {}
for line in open('/home/unix/csmillie/Gut_Human/csmillie/data/Sample_Info.txt'):
    line = line.rstrip().split()
    folder = line[0]
    sample = line[1]
    info[sample] = folder

# Extract info
[sample, site] = args.sid.split('.')
folder = info[sample]
site = sites[site]

# Get sequence IDs
seqs = {}
orgs = 'archaea bacteria fungi protozoa viral'.split()
for org in orgs:
    fn = './%s/%s_%s_bowtie2_contam.fastq' %(org, args.sid, org)
    for record in util.iter_fsq(fn):
        if len(record) > 0:
            seqs[record[0][1:]] = record
print 'Found %d sequences' %(len(seqs))

# Parse BAM file
out = open(args.out, 'w')
fn = '/home/unix/csmillie/Gut_Human/data/%s/%s/outs/possorted_genome_bam.bam' %(folder, site)
if not os.path.exists(fn):
    print 'BAM file not found'
for line in os.popen('samtools view -f 4 %s' %(fn)):
    read = line.split()[0]
    if read in seqs:

        cell = ''
        umi = ''

        if 'CB:Z' in line:
            cell = re.search('CB:Z:([ACGTacgtNn]*)', line).group(1)
        elif 'CR:Z' in line:
            cell = re.search('CR:Z:([ACGTacgtNn]*)', line).group(1)
        else:
            continue
        
        if 'UB:Z' in line:
            umi = re.search('UB:Z:([ACGTacgtNn]*)', line).group(1)
        elif 'UR:Z' in line:
            umi = re.search('UR:Z:([ACGTacgtNn]*)', line).group(1)
        else:
            continue
        
        record = seqs[read]
        record[0] = '%s;%s;%s' %(record[0], cell, umi)
        
        out.write('\n'.join(record) + '\n')

out.close()

# Remove UniVec sequences
cmd = 'sh 3.remove_univec.sh %s > %s.tmp; mv %s.tmp %s.uniclean' %(args.out, args.out, args.out, args.out)
os.system(cmd)
