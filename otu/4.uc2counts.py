import argparse

def parse_args():
    # Parse command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('-m', help='Dereplicated sequences mapping file', default='')
    parser.add_argument('-u', help='USEARCH output file', required=True)
    parser.add_argument('-s', help='Sample ID separator (regex)')
    parser.add_argument('-o', help='OTU table output file', required=True)
    args = parser.parse_args()
    return args


def map_derep(fn):
    # Map sequences to counts in each sample
    derep = {}
    for line in open(fn):
        sid, hits = line.rstrip().split('\t')
        derep[sid] = {}
        for hit in hits.split(' '):
            sa, count = hit.split(':')
            derep[sid][sa] = int(count)
    return derep


def map_uc(fn, derep={}, regex=''):
    # Convert USEARCH results to counts
    counts = {}
    samples = []
    for line in open(fn):
        if line.startswith('H'):
            line = line.rstrip().split()
            [sid, otu] = line[-2:]
            if otu not in counts:
                counts[otu] = {}
            if not derep:
                sa = re.search(regex, sid).group(1)
                if sa not in samples:
                    samples.append(sa)
                counts[otu][sa] = counts[otu].get(sa, 0) + 1
            else:
                for sa in derep[sid]:
                    if sa not in samples:
                        samples.append(sa)
                    counts[otu][sa] = counts[otu].get(sa, 0) + derep[sid][sa]
    return counts, samples


def write_otu_table(fn, counts, samples):
    
    # get samples and otus
    samples = sorted(samples)
    otus = sorted(counts.keys())
    
    # write otu matrix to file
    out = open(fn, 'w')
    out.write('sample_id\t' + '\t'.join(otus) + '\n')
    for sa in samples:
        outline = [sa]
        for otu in otus:
            try:
                outline.append('%d' %(counts[otu][sa]))
            except:
                outline.append('0')
        out.write('\t'.join(outline) + '\n')
    out.close()


# construct otu table from derep map and USEARCH results
args = parse_args()
if args.m:
    derep = map_derep(fn = args.m)
else:
    derep = {}
[counts, samples] = map_uc(fn=args.u, derep=derep, regex=args.s)
write_otu_table(args.o, counts, samples)
