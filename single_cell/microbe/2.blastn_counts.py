import argparse, glob, pdb, re, util
from Bio import pairwise2


# read input arguments
parser = argparse.ArgumentParser()
parser.add_argument('--fastq', help='FASTQ file', required=True)
parser.add_argument('--blast', help='BLAST file', required=True)
parser.add_argument('--sample', help='Sample name', required=True)
parser.add_argument('--prefix', help='OTU prefix', default='bacteria')
parser.add_argument('--out', help='Output prefix (counts)')
args = parser.parse_args()


# ---------------
# Initialize data
# ---------------
print('Initializing data')

# 1. Read FASTQ sequences into dictionary
seqs = {}
for record in util.iter_fsq(args.fastq):
    sid = record[0][1:]
    seq = record[1]
    seqs[sid] = seq

# 2. Map: microbial contigs to GCF IDs
contig2gcf = {}
for line in open('/home/unix/csmillie/aviv/db/refseq/meta/contig2gcf.txt'):
    contig, gcf = re.sub('"', '', line).rstrip().split('\t')
    contig2gcf[contig] = gcf        


# 3. Map: GCF IDs to taxonomy
gcf2sp = {}
gcf2gn = {}
for line in open('/home/unix/csmillie/aviv/db/refseq/meta/gcf.taxonomy_table.txt'):
    # header: (0=gcf, 1=name, 2=strain, 3=species id, 4=species, 5=genus id, 6=genus)
    if line.startswith('assembly_accession'):
        continue
    line = re.sub('"', '', line).rstrip().split('\t')
    gcf2sp[line[0]] = line[4]
    gcf2gn[line[0]] = line[6]


# --------------------------
# Check redundant reads/UMIs
# --------------------------


def check_umi(read, umi, seq):
    # This function checks whether a read/UMI/seq has already been counted
    # If the read or both the UMI and sequence have been seen before, it returns FALSE
    # Otherwise, it returns TRUE

    # Check if read has been mapped
    if read in umis:
        return False
    else:
        umis[read] = 1
    
    # Check if UMI/read has been mapped
    flag = True
    if umi in umis:
        for ref in umis[umi]:
            score = pairwise2.align.globalxx(seq, ref, score_only=True)
            if score >= .75*len(seq):
                flag = False
    else:
        umis[umi] = []
    
    # If not, then add to umis dict
    umis[umi].append(seq)
    return flag


# 1. Initialize variables
sp_count = {}
gn_count = {}
umis = {}
all_sp = {}
all_gn = {}
all_cells = {}

# 2. Parse BLAST report
print('Parsing BLAST report')
for line in open(args.blast):
    
    # Extract info
    line = line.rstrip().split('\t')
    read, cell, umi = line[0].split(';')
    umi = '%s.%s' %(cell, umi)
    acc = line[1]
    pct = float(line[2])
    seq = seqs[line[0]]
    gcf = contig2gcf[acc]
    evalue = float(line[10])
    sp = gcf2sp[gcf]
    gn = gcf2gn[gcf]
    
    # Append prefixes
    cell = '%s.%s' %(args.sample, cell)
    sp = '%s.%s' %(args.prefix, sp)
    gn = '%s.%s' %(args.prefix, gn)
    
    # Percent identity and E-value cutoff
    if pct <= 95 or evalue >= 1e-5:
        continue
    
    # Check if read or UMI has been mapped
    if not check_umi(read, umi, seq):
        continue
    
    # Count species and genus
    if cell not in sp_count:
        sp_count[cell] = {}
    if cell not in gn_count:
        gn_count[cell] = {}
    sp_count[cell][sp] = sp_count[cell].get(sp, 0) + 1
    gn_count[cell][gn] = gn_count[cell].get(gn, 0) + 1

    # Update rows/columns
    all_sp[sp] = 1
    all_gn[gn] = 1
    all_cells[cell] = 1

# Sort rows/columns
all_sp = sorted(all_sp)
all_gn = sorted(all_gn)
all_cells = sorted(all_cells)

print('Writing OTU tables')

# 3. Species table
out = open('%s.sp_counts.txt' %(args.out), 'w')
if len(all_cells) > 0:
    oi = ['OTU'] + all_cells
    out.write('\t'.join(oi) + '\n')
    for sp in all_sp:
        oi = [sp]
        for cell in all_cells:
            oi.append(sp_count[cell].get(sp, 0))
        out.write('\t'.join(map(str, oi)) + '\n')
out.close()

# 4. Genus table
out = open('%s.gn_counts.txt' %(args.out), 'w')
if len(all_cells) > 0:
    oi = ['OTU'] + all_cells
    out.write('\t'.join(oi) + '\n')
    for gn in all_gn:
        oi = [gn]
        for cell in all_cells:
            oi.append(gn_count[cell].get(gn, 0))
        out.write('\t'.join(map(str, oi)) + '\n')
out.close()
