import glob, re
from Bio import pairwise2

# Get all SAM files
sams = glob.glob('*/*.sam')

# Map contigs to GCFs
contig2gcf = {}
for line in open('/home/unix/csmillie/aviv/db/refseq/meta/contig2gcf.txt'):
    contig, gcf = line.rstrip().split('\t')
    contig2gcf[contig] = gcf

# Initialize data
counts = {}
umis = {}
gcf2reads = {}
samples = []
microbes = []

def test_umi(umi, seq):
    if umi not in umis:
        return True
    else:
        for ref in umis[umi]:
            score = pairwise2.align.globalxx(seq, ref)
            if score >= .75*len(seq):
                return False
        return True

# For every SAM file, get hits
for sam in sams:
    
    # Extract sample and host info
    [subject, site, host, suffix] = re.sub('.*\/', '', sam).split('.')
    
    # Parse mapped reads
    for line in open(sam):
        
        # Skip header
        if line.startswith('@'):
            continue
        
        # Get alignment info
        line = line.rstrip().split('\t')
        
        # Skip unmapped reads
        if line[1] == '4':
            continue

        # Bit flag is 0 or 16
        assert line[1] in ['0', '16']
        
        # Extract cell, umi, host, sequence
        [read, cell, umi] = line[0].split(';')
        acc = line[2]
        seq = line[9]
        gcf = contig2gcf[acc]

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
bad_gcfs = []

for i in range(len(all_gcfs)-1):
    gcf_i = all_gcfs[i]
    for j in range(i+1, len(all_gcfs)):
        gcf_j = all_gcfs[j]
        overlap = 1.*sum([1 for ri in gcf2reads[gcf_j] if ri in gcf2reads[gcf_i] else 0])/len(gcf2reads[gcf_j])
        if overlap >= .5:
            print '%s (size = %d) has %f overlap with %s (size = %d)' %(gcf_i, len(gcf2reads[gcf_i]), overlap, gcf_j, len(gcf2reads[gcf_j]))
            bad_gcfs.append(gcf_j)


print '\t' + '\t'.join(microbes)
for sample in samples:
    outline = [sample]
    for microbe in microbes:
        gcf = re.sub('.*?\\.', '', microbe)
        if gcf in bad_gcfs:
            continue
        outline.append(str(counts[sample].get(microbe, 0)))
    print '\t'.join(outline)
