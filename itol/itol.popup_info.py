import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--text2', help='input popup info (2 columns)')
parser.add_argument('--text3', help='input popup info (3 columns)')
args = parser.parse_args()

print '''
POPUP_INFO
SEPARATOR SPACE
DATA'''.strip()

if args.text2:
   for line in open(args.text2):
      line = line.rstrip().split()
      print line[0], 'info', line[1]

if args.text3:
   for line in open(args.text3):
      line = line.rstrip().split()
      print line[0], line[1], line[2]
