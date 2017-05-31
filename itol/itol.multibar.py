import argparse, random, tableau_colors

# get tableau colors
t20 = tableau_colors.colors['Tableau 20']
t20 = t20[0::2] + t20[1::2]

# parse command line arguments
parser = argparse.ArgumentParser()
parser.add_argument('-i', help='input tsv (1=leaf, 2=bar1, 3=bar2, etc)')
parser.add_argument('-p', help='convert negatives to positives', action='store_true', default=False)
parser.add_argument('-H', help='header', action='store_true', default=False)
parser.add_argument('-l', help='dataset label', default='label1')
args = parser.parse_args()

# get header info
h = open(args.i).readline().rstrip().split('\t')[1:]
if args.p:
   h = ['%s_neg' %(hi) for hi in h] + ['%s_pos' %(hi) for hi in h]
color = random.choice(t20)
field_colors = ' '.join(t20[:len(h)])
if args.H:
   field_labels = h
else:
   field_labels = map(str, range(1, len(h)+1))
field_labels = ' '.join(field_labels)

# print header
header = '''
DATASET_MULTIBAR
SEPARATOR SPACE
DATASET_LABEL %s
COLOR %s
FIELD_COLORS %s
FIELD_LABELS %s
DATA
''' %(args.l, color, field_colors, field_labels)
print header.strip()

# get barchart values
h = 0
for line in open(args.i):
   # skip header
   if args.H and h == 0:
      h = 1
      continue
   line = line.rstrip().split()
   if args.p == False:
      print ' '.join(line)
   else:
      out = [line[0]]
      for i in range(1, len(line)):
         if float(line[i]) < 0:
            out += [abs(float(line[i])), 0]
         else:
            out += [0, line[i]]
      print ' '.join(map(str, out))
