import argparse, os, re
import ssub3 as ssub

steps = 'demux cellranger cp dge tsne'.split()


# ---------------
# Input arguments
# ---------------

parser = argparse.ArgumentParser()
parser.add_argument('--folder', help='Illumina BCL run folder', default='')
parser.add_argument('--serial', help='Flowcell serial number', default='')
parser.add_argument('--fastq', help='FASTQ folder', default='')
parser.add_argument('--csv', help='CellRanger csv file (Lane, Sample, Index)', required=True)
parser.add_argument('--host', help='Host organism', default='hg19', choices=['hg19', 'mm10', 'hg19_and_mm10', 'mm10_premrna', 'hg19_premrna', 'hg19_and_mm10_premrna'])
parser.add_argument('--version', help='CellRanger version', default='3.0.2', choices=['1.1.0', '1.3.1', '2.1.1', '3.0.2'])
parser.add_argument('--enter' ,help='Enter pipeline at step', default='demux', choices=steps)
parser.add_argument('--exit', help='Exit pipeline at step', default='tsne', choices=steps)
parser.add_argument('-m', help='Memory (default=32)', type=int, default=32)
parser.add_argument('-p', help='Print out commands?', default=False, action='store_true')
parser.add_argument('-w', help='Write scripts and exit', default=False, action='store_true')
args = parser.parse_args()


# --------------------
# Initialize variables
# --------------------

# Make output folder
if not os.path.exists('outs'):
    os.mkdir('outs')

ref = '/home/unix/csmillie/Gut_Human/csmillie/gtf/%s' %(args.host)
sigdb = '/home/unix/csmillie/db/immgen.%s.db.txt' %(args.host)

# Get flowcell from RunInfo file
try:
    run_info = os.path.join(args.folder, 'RunInfo.xml')
    run_info = '\n'.join(open(run_info).readlines())
    serial = re.search('<Flowcell>(.*?)</Flowcell>', run_info).group(1)
except:
    serial = args.serial
if not serial:
    quit('No flowcell serial')

# Extract sample index from filename
print('Predicted serial: %s' %(serial))
if not args.fastq:
    fastq_path = './%s/outs/fastq_path/' %(serial)
else:
    fastq_path = args.fastq
print('Predicted FASTQ path: %s' %(fastq_path))

# Get sample names
samples = [line.rstrip().split(',')[1] for line in open(args.csv) if not line.startswith('Lane')]

def check_step(name):
    if steps.index(args.enter) <= steps.index(name) and steps.index(args.exit) >= steps.index(name):
        return True

pipeline = []


# -----------
# Demultiplex
# -----------

demux_fn = '/home/unix/csmillie/code/10X/wrap_10x_demux.sh'
if not os.path.exists(demux_fn):
    quit('Error: could not find %s' %(demux_fn))
cmd1 = "cat %s | sed -e 's INDIR %s ' -e 's CSV %s ' -e 's VERSION %s ' > ./demux.sh" %(demux_fn, args.folder, args.csv, args.version)
os.system(cmd1)
if check_step('demux'):
    pipeline.append(['sh ./demux.sh'])


# -----------
# Cell Ranger
# -----------

run10x_fn = '/home/unix/csmillie/code/10X/wrap_10x_run_parallel.sh'
cmd2 = "cat %s | sed -e 's CSV %s ' -e 's REF %s ' -e 's SERIAL %s ' -e 's VERSION %s ' -e 's NUM %s ' -e 's FASTQ_PATH %s ' > ./run10x.sh" %(run10x_fn, args.csv, ref, serial, args.version, len(samples), fastq_path)
if not os.path.exists(run10x_fn):
    quit('Error: could not find %s' %(run10x_fn))
os.system(cmd2)
if check_step('cellranger'):
    pipeline.append(['sh ./run10x.sh'])

if args.w:
    quit()

# -----------------
# Copy output files
# -----------------

if check_step('cp'):
    cmds = []
    if not os.path.exists('./outs/CellRanger'):
        os.mkdir('./outs/CellRanger')
    for sample in samples:    
        # Make sample directory
        if not os.path.exists('./outs/CellRanger/%s' %(sample)):
            os.mkdir('./outs/CellRanger/%s' %(sample))    
        # Copy folders
        cmdi = 'cp -r ./%s/outs/*gene_bc_matrices ./outs/CellRanger/%s/' %(sample, sample)
        cmdk = 'cp ./%s/outs/web_summary.html ./outs/CellRanger/%s/' %(sample, sample)
        cmdl = 'cp ./%s/outs/metrics_summary.csv ./outs/CellRanger/%s/' %(sample, sample)
        cmds.append('%s; %s; %s' %(cmdi, cmdk, cmdl))
    pipeline.append(cmds)
    pipeline.append(['tar -cvzf ./outs/CellRanger.tgz --remove-files -C ./outs ./CellRanger'])


# --------------
# Construct DGEs
# --------------

if check_step('dge'):
    
    # DGE for individual samples
    cmds = []
    if not os.path.exists('./outs/DGE'):
        os.mkdir('./outs/DGE')
    samples = [line.rstrip().split(',')[1] for line in open(args.csv) if not line.startswith('Lane')]
    for sample in samples:
        folder = './%s/outs/filtered_gene_bc_matrices/%s' %(sample, args.host)
        cmdi = 'Rscript /home/unix/csmillie/code/10X/make_counts.r --folder %s --prefix %s --out ./outs/DGE/%s.dge.txt' %(folder, sample, sample)
        cmdj = 'gzip ./outs/DGE/%s.dge.txt' %(sample)
        cmds.append('%s; %s' %(cmdi, cmdj))
    pipeline.append(cmds)
    
    # DGE for all samples
    out = open('all.map.txt', 'w')
    for sample in samples:
        out.write('NA\t./outs/DGE/%s.dge.txt.gz\t.*\n' %(sample))
    out.close()
    pipeline.append(['Rscript /home/unix/csmillie/code/single_cell/merge_dges.r --map all.map.txt --ming 100 --out ./outs/DGE/all.dge.txt; gzip ./outs/DGE/all.dge.txt'])
    

# --------------
# Calculate TSNE
# --------------

if check_step('tsne'):
    cmds = []
    if not os.path.exists('./outs/Portal'):
        os.mkdir('./outs/Portal')
    for sample in samples:
        cmds.append('Rscript /home/unix/csmillie/code/10X/do_seurat.r --name ./outs/Portal/%s --dge ./outs/DGE/%s.dge.txt.gz --ming 500 --sigdb %s' %(sample, sample, sigdb))
    cmds.append('Rscript /home/unix/csmillie/code/10X/do_seurat.r --name ./outs/Portal/all --dge ./outs/DGE/all.dge.txt.gz --ming 500 --sigdb %s --out' %(sigdb))
    pipeline.append(cmds)


# ------------
# Run pipeline
# ------------

import ssub3 as ssub
Ssub = ssub.Submitter()
Ssub.m = args.m
ssub.p = args.p

for group in pipeline:
    print('\n'.join(group))

for cmd in pipeline:
    Ssub.add_task(cmd)

Ssub.submit()
