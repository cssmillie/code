import argparse, os, re, util

# get UNMAPPED reads from bam file (-F for mapped reads)

# get input arguments
parser = argparse.ArgumentParser()
parser.add_argument('--bam', help='Input BAM file (.bam)', required=True)
parser.add_argument('--out', help='Output FASTQ file (.fastq)', required=True)
parser.add_argument('--subset', help='Sequence IDs to subset (.txt)', default='')
parser.add_argument('--add_barcodes', help='Add cell/UMI barcodes?', default=False, action='store_true')
args = parser.parse_args()

# check arguments
if not args.subset and not args.add_barcodes:
    cmd = 'samtools bam2fq -f 4 %s > %s' %(args.bam, args.out)
    print(cmd)
    quit()

# subset sequence IDs
if args.subset:
    subset = {li.rstrip(): 1 for li in open(args.subset)}
else:
    subset = {}

# parse BAM file
out = open(args.out, 'w')
for line in os.popen('samtools view -f 4 %s' %(args.bam)):

    # parse sam record
    read = line.split()[0]
    cell = umi = ''

    # subset reads
    if args.subset:
        if read not in subset:
            continue

    # add barcodes
    if args.add_barcodes:
    
        # add cell barcode
        if 'CB:Z' in line:
            cell = re.search('CB:Z:([ACGTacgtNn]*)', line).group(1)
        elif 'CR:Z' in line:
            cell = re.search('CR:Z:([ACGTacgtNn]*)', line).group(1)
        else:
            continue

        # add umi barcode
        if 'UB:Z' in line:
            umi = re.search('UB:Z:([ACGTacgtNn]*)', line).group(1)
        elif 'UR:Z' in line:
            umi = re.search('UR:Z:([ACGTacgtNn]*)', line).group(1)
        elif 'RX:Z' in line:
            umi = re.search('RX:Z:([ACGTacgtNn]*)', line).group(1)
        else:
            continue
        
        # get sequence info
        sid = '%s;%s;%s' %(read, cell, umi)

    else:
        sid = read

    # extract line info
    line = line.rstrip().split()
    seq = line[9]
    qua = line[10]
    
    # write FASTQ record
    out.write('@%s\n%s\n+\n%s\n' %(sid, seq, qua))

out.close()
