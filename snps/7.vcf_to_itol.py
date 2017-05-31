import argparse

# get input arguments
parser = argparse.ArgumentParser()
parser.add_argument('-i', help='input vcf file')
parser.add_argument('-m', help='number of total alleles per sample', default=0, type=int)
parser.add_argument('-M', help='number of minor alleles per sample', default=0, type=int)
parser.add_argument('-n', help='number of total alleles per snp', default=0, type=int)
parser.add_argument('-N', help='number of minor alleles per snp', default=0, type=int)
parser.add_argument('-o', help='output prefix')
parser.add_argument('-c', help='map of samples to colors')
args = parser.parse_args()

# parse vcf
cmd = 'python ~/aviv/box/snps/parse_vcf.py --vcf %s --out %s.dna' %(args.i, args.i)
print cmd

# filter dna matrix
cmd = 'Rscript ~/aviv/box/snps/filter_dna.r --in %s.dna --m %d --M %d --n %d --N %d --out %s.f.dna' %(args.i, args.m, args.M, args.n, args.N, args.o)
print cmd

# convert to fasta
cmd = 'cat %s.f.dna | sed "s/^/>/g" | tr ";" "\\n" > %s.fasta' %(args.o, args.o)
print cmd

# make phylogeny
cmd = 'FastTree -nt -gtr < %s.fasta > %s.newick' %(args.o, args.o)
print cmd

# make itol plots
cmd = 'python ~/box/itol.color_strips.py -i %s --label cell_types' %(args.c)
print cmd


