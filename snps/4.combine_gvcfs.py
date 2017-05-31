import os.path, ssub

ref = '/broad/mccarroll/software/metadata/merged_reference/hg19_mm10_transgenes/hg19_mm10_transgenes.fasta'

submitter = ssub.Submitter()

# Map intervals -> GVCF files
def get_intervals():
    intervals = {}
    for i in range(100):
        fns = []
        collections = [line.rstrip() for line in open('collections.txt')]
        for collection in collections:
            barcodes = [line.rstrip() for line in open('%s.barcodes.txt' %(collection))]
            for barcode in barcodes:
                fns.append('./vcf/%s/%s/%s_%s.%d.g.vcf' %(collection, barcode, collection, barcode, i))
        intervals[i] = fns
    return intervals

intervals = get_intervals()

# Hierarchically combine GVCFs
all = {}
iteration = 0
while True:
    iteration += 1
    # Stop when there are no more files to merge
    if sum([len(fns) for fns in intervals.values()]) == len(intervals):
        break
    # After each iteration, we submit these commands to the queue
    cmds = []
    # Get merge commands for every interval
    for i in intervals:
        fns = intervals[i]
        count = 0
        if len(fns) <= 1:
            continue
        else:
            new_fns = []
            while len(fns) > 0:
                count += 1
                j = min(200, len(fns))
                v = ' '.join(['-V %s' %(fn) for fn in fns[:j]])
                fns = fns[j:]
                new_fn = './vcf/combine/interval%d.iter%d.%d.combine.g.vcf' %(i, iteration, count)
                new_fns.append(new_fn)
                all[new_fn] = 1
                cmd = 'java -jar ~/bin/GenomeAnalysisTK.jar -T CombineGVCFs -R %s %s -o %s' %(ref, v, new_fn)
                cmds.append(cmd)
            intervals[i] = new_fns
    while len(cmds) > 0:
        # Only run commands that haven't finished
        new_cmds = []
        for cmd in cmds:
            new_fn = cmd.split()[-1]
            idx_fn = new_fn + '.idx'
            if not os.path.exists(idx_fn):
                new_cmds.append(cmd)
        cmds = new_cmds
        if len(cmds) > 0:
            submitter.submit_and_wait(cmds)

# Get files to keep/remove
u = []
v = all.keys()
for i in intervals:
    fns = intervals[i]
    if len(fns) == 1:
        u.append(fns[0])
        v.remove(fns[0])
    else:
        print intervals
        quit('error, length > 1')

# Write these to file
out = open('./vcf/combine/keep.txt', 'w')
for ui in u:
    out.write(ui + '\n')
out.close()

out = open('./vcf/combine/remove.txt', 'w')
for vi in v:
    out.write(vi + '\n')
out.close()
        
