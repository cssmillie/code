import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--colors', help='input color map (leaf -> hex)')
parser.add_argument('--groups', help='input group map (leaf -> group)')
args = parser.parse_args()

print '''
DATASET_COLORSTRIP
SEPARATOR SPACE
DATASET_LABEL label1
COLOR #000000
DATA'''.strip()

if args.colors:
   for line in open(args.colors):
      line = line.rstrip().split()
      print line[0], line[1]

if args.groups:
   import tableau_colors
   t20 = tableau_colors.colors['Tableau 20']
   t20 = t20[0::2] + t20[1::2]
   g2c = {}
   for line in open(args.groups):
      line = line.rstrip().split()
      if line[1] not in g2c:
         g2c[line[1]] = t20[0]
         t20 = t20 + [t20.pop(0)]
      print line[0], g2c[line[1]]
