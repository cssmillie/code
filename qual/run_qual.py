import argparse, os.path, ssub

# Initialize submitter
Submitter = ssub.Submitter()

# Get input arguments
parser = argparse.ArgumentParser()
parser.add_argument('--folder', help='folder name')
parser.add_argument('--list', help='folder list')
parser.add_argument('--field', help='collection field')
parser.add_argument('--out', help='output prefix')
args = parser.parse_args()

# Get folders
if args.list:
    folders = [line.rstrip() for line in open(args.list)]
else:
    folders = [args.folder]
print '\n\nFound folders:\n' + '\n'.join(folders)

# Get files
dsq_folders = [folder for folder in folders if os.path.exists('/home/unix/csmillie/aviv/data/%s/raw/QC_files' %(folder))]
print '\n\nDSQ folders: ' + '\n'.join(dsq_folders)
dge_folders = [folder for folder in folders if os.path.exists('/home/unix/csmillie/aviv/data/%s/dge' %(folder))]
print '\n\nDGE folders: ' + '\n'.join(dge_folders)

# Run jobs
dsq_cmds = ['Rscript /home/unix/csmillie/aviv/code/qual/dsq_qual.r --folder %s --out %s' %(folder, folder) for folder in dsq_folders]
dge_cmds = ['Rscript /home/unix/csmillie/aviv/code/qual/dge_qual.r --folder %s --field 1 --out %s' %(folder, folder) for folder in dge_folders]
c_dsq_outs = ['%s.collections.dsq_qual.txt' %(folder) for folder in dsq_folders]
e_dsq_outs = ['%s.experiments.dsq_qual.txt' %(folder) for folder in dsq_folders]
c_dge_outs = ['%s.collections.dge_qual.txt' %(folder) for folder in dge_folders]
e_dge_outs = ['%s.experiments.dge_qual.txt' %(folder) for folder in dge_folders]

# Submit
cmds = dsq_cmds + dge_cmds
outs = c_dsq_outs + e_dsq_outs + c_dge_outs + e_dge_outs
print 'Running %d jobs with 1 retry' %(len(cmds))
print 'Example job: %s' %(cmds[0])
print 'Example out: %s' %(outs[0])
Submitter.submit_and_wait(cmds, outs=outs, retry=1)

# Function to concatenate files
def concat_fns(fns, out_fn):
    header = open(fns[0]).readline()
    out = open(out_fn, 'w')
    out.write(header)
    for fn in fns:
        h = 0
        for line in open(fn):
            if h == 0:
                h += 1
                continue
            out.write(line)
    out.close()

# Concatenate dsq outs
print 'Concatenating output files'
concat_fns(c_dsq_outs, '%s.collections.dsq_qual.txt' %(args.out))
concat_fns(e_dsq_outs, '%s.experiments.dsq_qual.txt' %(args.out))
concat_fns(c_dge_outs, '%s.collections.dge_qual.txt' %(args.out))
concat_fns(e_dge_outs, '%s.experiments.dge_qual.txt' %(args.out))
