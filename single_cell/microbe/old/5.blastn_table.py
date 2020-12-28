import glob, pdb, re, util
from Bio import pairwise2

# Get all BLASTN files
sams = glob.glob('*/*.blastn')

# Get all sequences
seqs = {}
for fn in glob.glob('*.contam.fastq'):
    for record in util.iter_fsq(fn):
        seqs[record[0][1:]] = record[1]

# Initialize data
counts = {}
umis = {}
samples = []
microbes = []
gcf2reads = {}

# Map contigs to GCFs
contig2gcf = {}
for line in open('/home/unix/csmillie/aviv/db/refseq/meta/contig2gcf.txt'):
    contig, gcf = line.rstrip().split('\t')
    contig2gcf[contig] = gcf        

# Map GCF to genome name
gcf2name = {}
for line in open('/home/unix/csmillie/aviv/db/refseq/meta/acc2name.txt'):
    if line.startswith('#'):
        continue
    line = line.rstrip().split('\t')
    gcf2name[line[0]] = line[1]
    
def test_umi(umi, seq):
    if umi not in umis:
        return True
    else:
        for ref in umis[umi]:
            score = pairwise2.align.globalxx(seq, ref)
            if score >= .75*len(seq):
                return False
        return True

# Get best hits
blast = {}
uid = 0
for sam in sams:
    [subject, site, host, suffix] = re.sub('.*\/', '', sam).split('.')
    for line in open(sam):
        line = line.rstrip().split('\t')
        if float(line[2]) <= 95:
            continue
        line.append(sam)
        read = line[0]
        evalue = float(line[10])
        blast[uid] = line
        uid += 1


# For every SAM file, get hits
print 'Parsing BLAST results'
for uid in blast:
    
    line = blast[uid]
    
    # Extract sample and host info
    [subject, site, host, suffix] = re.sub('.*\/', '', line[-1]).split('.')
    
    # Percent identity & E-value cutoff
    if float(line[10]) >= 1e-5:
        continue
        
    # Extract cell, umi, host, sequence
    [read, cell, umi] = line[0].split(';')
    acc = line[1]
    seq = seqs[line[0]]

    try:
        gcf = contig2gcf[acc]
    except:
        continue
    
    # Keep track of reads mapped to each GCF
    if gcf not in gcf2reads:
        gcf2reads[gcf] = {}
        gcf2reads[gcf][read] = 1
    
    # Test if UMI has already been mapped
    umi = '%s.%s' %(cell, umi)
    if test_umi(umi, seq):
        
        sample = '%s.%s.%s' %(subject, site, cell)
        microbe = '%s.%s' %(host, gcf)
        
        if sample not in samples:
            samples.append(sample)            
        if microbe not in microbes:
            microbes.append(microbe)

        if sample not in counts:
            counts[sample] = {}
        if microbe not in counts[sample]:
            counts[sample][microbe] = 0
        counts[sample][microbe] += 1


# Remove redundant GCFs
all_gcfs = sorted(gcf2reads.keys(), key=lambda x: len(gcf2reads[x]), reverse=True)
bad_gcfs = {}

print 'Removing redundant GCFs'
for i in range(len(all_gcfs)-1):
    gcf_i = all_gcfs[i]
    if gcf_i in bad_gcfs:
        continue
    for j in range(i+1, len(all_gcfs)):
        gcf_j = all_gcfs[j]
        if gcf_j in bad_gcfs:
            continue
        overlap = 1.*sum([1 if ri in gcf2reads[gcf_i] else 0 for ri in gcf2reads[gcf_j]])/len(gcf2reads[gcf_j])
        if overlap >= .5:
            print '%s (size = %d) has %f overlap with %s (size = %d)' %(gcf_i, len(gcf2reads[gcf_i]), overlap, gcf_j, len(gcf2reads[gcf_j]))
            print 'combining %s <-- %s' %(gcf2name[gcf_i], gcf2name[gcf_j])
            bad_gcfs[gcf_j] = 1


old_size = len(microbes)
microbes = [microbe for microbe in microbes if re.sub('^.*?\.', '', microbe) not in bad_gcfs]
new_size = len(microbes)
print 'removed redundant %d microbes' %(old_size - new_size)

out = open('test', 'w')
out.write('\t' + '\t'.join(microbes) + '\n')
for sample in samples:
    outline = [sample]
    for microbe in microbes:
        outline.append(str(counts[sample].get(microbe, 0)))
    out.write('\t'.join(outline) + '\n')
out.close()
