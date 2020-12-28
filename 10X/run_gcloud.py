
import argparse, os, re

# parse input arguments
parser = argparse.ArgumentParser()
parser.add_argument('--folder', help='demux folder')
parser.add_argument('--csv', help='input samples csv')
args = parser.parse_args()

def get_cpath(fn):
    # get path on google cloud
    return '%s%s' %(os.environ['CLOUD'], os.path.abspath(fn))

# copy file to google cloud
folder = os.path.abspath(args.folder)
cmd = 'gsutil -m rsync -r %s %s' %(folder, get_cpath(folder))

# samples.csv file for cloud
out = open('samples.gcloud.csv', 'w')
out.write('Sample,Reference,Flowcell,Lane,Index,Chemistry,DataType,FeatureBarcodeFile')
for line in open(args.csv):
    if line.startswith('Lane'):
        continue
    lane, sample, index = line.rstrip().split(',')
    



h = 'Sample,Reference,Flowcell,Lane,Index,Chemistry,DataType,FeatureBarcodeFile'
for line in open(args.csv):

