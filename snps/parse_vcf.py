import argparse, random, re, scipy.stats
from collections import Counter

parser = argparse.ArgumentParser()
parser.add_argument('--vcf', help='input vcf', required=True)
parser.add_argument('--out', help='output dna matrix', required=True)
args = parser.parse_args()

header = []

out = open(args.out, 'w')

for line in open(args.vcf):

    if line.startswith('#CHROM'):
        # print header
        header = '\t'.join(line.rstrip().split('\t')[9:])
        out.write(header + '\n')
    
    if line.startswith('HUMAN'):
        line = line.rstrip().split('\t')
        
        # get allele counts
        ac = map(float, re.search('AC=(.*?);', line[7]).group(1).split(','))
        tc = float(re.search('AN=(.*?);', line[7]).group(1))
        allele_counts = [tc-sum(ac)] + ac
        
        # get major allele
        temp = sorted([(v,i) for i,v in enumerate(allele_counts)])
        Mc, Mi = temp[-1]
        mc, mi = temp[-2]
        
        # test homozygous
        hom = 1 - scipy.stats.binom.sf(Mc-1, tc, .97)
        # test heterozygous 
        het = scipy.stats.binom_test(Mc, tc, .50)

        # get expected alleles
        if hom > .01 and het > .01:
            continue
        # homozygous
        if hom > .01 and het < .01:
            ok = [Mi]
        # heterozygous
        elif hom < .01 and het > .01:
            ok = [Mi, mi]
        # skewed freqs
        else:
            ok = [Mi, mi]

        # randomize data
        nts = 'A C G T'.split()
        random.shuffle(nts)
        
        # print genotypes
        dna = []
        for i in range(9, len(line)):
            entry = line[i]
            # write missing data
            if entry.startswith('./.'):
                dna.append('-')
                continue
            # call cell genotype
            alleles = map(int, entry.split(':')[0].split('/'))
            di = nts[0]
            for allele in alleles:
                if allele not in ok:
                    di = nts[allele]
            dna.append(di)
        # skip monomorphic sites
        if len(set([di for di in dna if di != '-'])) <= 1:
            continue
        # write out
        out.write('\t'.join(dna) + '\n')
out.close()
