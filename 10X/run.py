import argparse, os, re, ssub

steps = 'demux cellranger cp dge tsne'.split()

# Get input arguments
parser = argparse.ArgumentParser()
parser.add_argument('--folder', help='data directory', required=True)
parser.add_argument('--barcodes', help='sample -> barcodes (tsv)', required=True)
parser.add_argument('--human', help='human?', action='store_true', default=False)
parser.add_argument('--mouse', help='mouse?', action='store_true', default=False)
parser.add_argument('--version', help='cellranger version', default='1.3.1', choices=['1.1.0', '1.3.1'])
parser.add_argument('--enter' ,help='enter pipeline at step', default='demux', choices=steps)
parser.add_argument('--exit', help='exit pipeline at step', default='tsne', choices=steps)
parser.add_argument('-m', help='Memory (default=32)', type=int, default=32)
parser.add_argument('-p', help='print out commands?', default=False, action='store_true')
args = parser.parse_args()

# Make output folder
if not os.path.exists('outs'):
    os.mkdir('outs')

# Get reference genome
if args.human:
    host = 'hg19'
elif args.mouse:
    host = 'mm10'
else:
    quit('Must specify human or mouse')
ref = '/seq/regev_genome_portal/SOFTWARE/10X/refdata-cellranger-1.2.0/refdata-cellranger-%s-1.2.0' %(host)
sigdb = '/home/unix/csmillie/db/immgen.%s.db.txt' %(host)

# Get flowcell from RunInfo file
run_info = os.path.join(args.folder, 'RunInfo.xml')
run_info = '\n'.join(open(run_info).readlines())
serial = re.search('<Flowcell>(.*?)</Flowcell>', run_info).group(1)

# Extract sample index from filename
print 'Predicted serial: %s' %(serial)
fastq_path = './%s/outs/fastq_path/' %(serial)
print 'Predicted FASTQ path: %s' %(fastq_path)

# Get sample names
samples = [line.rstrip().split()[0] for line in open(args.barcodes)]

def check_step(name):
    if steps.index(args.enter) <= steps.index(name) and steps.index(args.exit) >= steps.index(name):
        return True

pipeline = []

# Demultiplex
if not os.path.exists('/home/unix/csmillie/code/10X/wrap_10x_demux.sh'):
    quit('Error: could not find /home/unix/csmillie/code/10X/wrap_10x_demux.sh')
cmd1 = "cat /home/unix/csmillie/code/10X/wrap_10x_demux.sh | sed -e 's INDIR %s ' -e 's VERSION %s ' > ./demux.sh" %(args.folder, args.version)
os.system(cmd1)
if check_step('demux'):
    pipeline.append(['sh ./demux.sh'])

# Cell Ranger
cmd2 = "cat /home/unix/csmillie/code/10X/wrap_10x_run_parallel.sh | sed -e 's BARCODES %s ' -e 's REF %s ' -e 's FQPATH %s ' -e 's VERSION %s ' > ./run10x.sh" %(args.barcodes, ref, fastq_path, args.version)
if not os.path.exists('/home/unix/csmillie/code/10X/wrap_10x_run_parallel.sh'):
    quit('Error: could not find /home/unix/csmillie/code/10X/wrap_10x_run_parallel.sh')
os.system(cmd2)
if check_step('cellranger'):
    pipeline.append(['sh ./run10x.sh'])

# Copy Cell Ranger output
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
        cmdj = 'cp -r ./%s/outs/analysis ./outs/CellRanger/%s/' %(sample, sample)
        cmdk = 'cp ./%s/outs/web_summary.html ./outs/CellRanger/%s/' %(sample, sample)
        cmdl = 'cp ./%s/outs/metrics_summary.csv ./outs/CellRanger/%s/' %(sample, sample)
        cmds.append('%s; %s; %s; %s' %(cmdi, cmdj, cmdk, cmdl))
    pipeline.append(cmds)
    pipeline.append(['tar -cvzf ./outs/CellRanger.tgz --remove-files -C ./outs ./CellRanger'])
    

# Construct DGEs
if check_step('dge'):
    
    # DGE for individual samples
    cmds = []
    if not os.path.exists('./outs/DGE'):
        os.mkdir('./outs/DGE')
    samples = [line.rstrip().split()[0] for line in open(args.barcodes)]
    for sample in samples:
        folder = './%s/outs/filtered_gene_bc_matrices/%s' %(sample, host)
        cmdi = 'Rscript /home/unix/csmillie/code/10X/make_counts.r %s %s ./outs/DGE/%s.dge.txt' %(folder, sample, sample)
        cmdj = 'gzip ./outs/DGE/%s.dge.txt' %(sample)
        cmds.append('%s; %s' %(cmdi, cmdj))
    pipeline.append(cmds)
    
    # DGE for all samples
    out = open('all.map.txt', 'w')
    for sample in samples:
        out.write('NA\t./outs/DGE/%s.dge.txt.gz\t.*\n' %(sample))
    out.close()
    pipeline.append(['Rscript /home/unix/csmillie/code/merge_dges.r --map all.map.txt --ming 500 --out ./outs/DGE/all.dge.txt; gzip ./outs/DGE/all.dge.txt'])
    
# Calculate TSNE
if check_step('tsne'):
    cmds = []
    if not os.path.exists('./outs/Portal'):
        os.mkdir('./outs/Portal')
    for sample in samples:
        cmds.append('Rscript /home/unix/csmillie/code/10X/do_seurat.r --name ./outs/Portal/%s --dge ./outs/DGE/%s.dge.txt.gz --ming 500 --sigdb %s' %(sample, sample, sigdb))
    cmds.append('Rscript /home/unix/csmillie/code/10X/do_seurat.r --name ./outs/Portal/all --dge ./outs/DGE/all.dge.txt.gz --ming 500 --sigdb %s --out' %(sigdb))
    pipeline.append(cmds)

# Run pipeline
for cmd in pipeline:
    print '\n'.join(cmd)

import ssub
Ssub = ssub.Submitter()
Ssub.m = args.m
ssub.p = args.p
Ssub.q = 'long'
Ssub.submit_pipeline(pipeline)
