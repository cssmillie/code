import argparse, os, re, ssub3, sys


# ---------
# filenames
# ---------
#
# Unmapped, polyA-trimmed, complexity filtered, non-human reads
# DEL ./*fastq = unmapped reads from bam file
# DEL ./*.good.tmp.fastq = complexity-filtered reads
# DEL ./*.trim.fastq = trimmed adapters and min length > 20
# DEL ./nonhuman/*.fastq = non-human reads
# DEL ./nonhuman/*.human_bowtie2_contam.fastq = human reads
# DEL ./nonhuman/*.log = kneaddata logfile
# DEL ./nonhuman/*.univec = univec sequences
# ./nonhuman/*.clean.fastq = trimmed, filtered, non-human, non-univec sequences
#
# Kraken output
# ./kraken/*.kraken.report.txt
# ./kraken/*.kraken.out
# ./kraken/*.bracken
#
# Kneaddata output
# ./org/*.org.fastq = organism fastq files
# ./org/*.org.blastn = blastn output


# ---------------
# input arguments
# ---------------

# parse arguments
parser = argparse.ArgumentParser()
parser.add_argument('--fastq', help='FASTQ file', default='')
parser.add_argument('--bam', help='BAM file', default='')
parser.add_argument('--prefix', help='output prefix', default='')
parser.add_argument('--method', help='method', choices=['knead', 'kraken', 'both'], default='both')
parser.add_argument('--kraken_db', help='kraken database', default='/home/unix/csmillie/aviv/db/kraken')
parser.add_argument('--step', help='pipeline step 1-10 (see code comments)', default='0')
parser.add_argument('--overwrite', help='overwrite files?', default=False, action='store_true')
parser.add_argument('-m', help='memory', type=int, default=64)
parser.add_argument('-t', help='runtime', default='8:00:00')
parser.add_argument('-p', help='print', default=False, action='store_true')
args = parser.parse_args()

# fix arguments
if args.fastq and args.bam:
    quit('error: must provide FASTQ or BAM')
if not args.fastq and not args.bam:
    quit('error: must provide FASTQ or BAM')
if not args.prefix:
    if args.fastq:
        args.prefix = os.path.basename(args.fastq)
        args.prefix = os.path.splitext(args.prefix)[0]
    else:
        args.prefix = os.path.basename(args.bam)
        args.prefix = os.path.splitext(args.prefix)[0]

# step range
args.step = args.step.split(',')
if len(args.step) == 1:
    args.step = [int(args.step[0]), float('inf')]
else:
    args.step = [int(args.step[0]), int(args.step[1])]

# print arguments
if not args.p:
    print('FASTQ: %s\nBAM: %s\nPrefix: %s\nSteps: %s' %(args.fastq, args.bam, args.prefix, args.step))


# --------------
# job management
# --------------

Submitter = ssub3.Submitter()
Submitter.m = 64
Submitter.t = '12:00:00'
Submitter.o = 'microbe'

def check_step(step, outs=[]):
    if type(outs) == type(''):
        outs = [outs]
    if step < args.step[0] or step > args.step[1]:
        if not args.p:
            sys.stderr.write('Skipping step %s\n' %(step))
        return False
    if len(outs) == 0 or args.overwrite:
        if not args.p:
            sys.stderr.write('Running step %s\n' %(step))
        return True
    for out in outs:
        if os.path.exists(out) and os.path.getsize(out) > 0:
            pass
        else:
            if not args.p:
                sys.stderr.write('"%s" does not exist\n' %(out))
            break
    else:
        if not args.p:
            sys.stderr.write('Skipping step %s\n' %(step))
        return False
    if not args.p:
        sys.stderr.write('Running step %s\n' %(step))
    return True


# ------------------------------------
# pipeline step 1: get non-human reads
# ------------------------------------

cmds = ''

# step 1. get unmapped reads
# - input: bam file
# - output: fastq file (*.fastq)
if check_step(1, './nonhuman/%s.clean.fastq' %(args.prefix)):
    if args.bam:
        cmds += 'python /home/unix/csmillie/code/seq/bam_to_fastq.py --bam %s --out %s.fastq --add_barcodes; ' %(args.bam, args.prefix)
        args.fastq = '%s.fastq' %(args.prefix)

# step 2. remove polyA tails and filter low complexity sequences
# - input: %s.fastq (fastq file)
# - output:
#   - %s.good.tmp.fastq (trimmed and filtered sequences)
#   - %s.bad.tmp.fastq (low complexity sequences)
if check_step(2, './nonhuman/%s.clean.fastq' %(args.prefix)):
    cmds += 'prinseq-lite -fastq %s -lc_method dust -lc_threshold 16 -trim_tail_right 6 -trim_tail_left 6 -out_good %s.good.tmp -out_bad %s.bad.tmp; ' %(args.fastq, args.prefix, args.prefix)

# step 2. trim adapter sequences and require read length >= 20
if check_step(2, './nonhuman/%s.clean.fastq' %(args.prefix)):
    cmds += 'cutadapt -g "file:/home/unix/csmillie/bin/cutadapt.5p.TruSeq3-PE-2-SS2.fa" -m 20 -o %s.trim.fastq %s.good.tmp.fastq; ' %(args.prefix, args.prefix)
    
# step 3. cleanup files to save space
# - remove:
#   - %s.fastq
#   - %s.good.tmp.fastq, %s.bad.tmp.fastq
# - keep:
#   - %s.trim.fastq (trimmed and length filtered sequences)
if check_step(3, './nonhuman/%s.clean.fastq' %(args.prefix)):
    cmds += 'rm %s.fastq' %(args.prefix)
    cmds += 'rm %s.good.tmp.fastq; ' %(args.prefix)
    cmds += 'rm %s.bad.tmp.fastq; ' %(args.prefix)

# step 4. remove human reads
# - input: %s.trim.fastq (trimmed and filtered sequences)
# - output:
#   - ./nonhuman/%.fastq: non-human reads
#   - ./nonhuman/%s_*human_bowtie2_contam.fastq: human reads
#   - ./nonhuman/%s.log: logfile
if check_step(4, './nonhuman/%s.clean.fastq' %(args.prefix)):
    cmds += 'kneaddata --input %s.trim.fastq -db /home/unix/csmillie/aviv/db/bowtie2/human/human --bypass-trim --output nonhuman --output-prefix %s; ' %(args.prefix, args.prefix)

# step 5. cleanup files to save space
# - remove:
#   - %s.trim.fastq: trimmed non-mapped reads
#   - ./nonhuman/%s_*human_bowtie2_contam.fastq: human reads
#   - ./nonhuman/%s.log: logfile
# - keep:
#   - ./nonhuman/%s.fastq: non-human reads
if check_step(5, './nonhuman/%s.clean.fastq' %(args.prefix)):
    cmds += 'rm %s.trim.fastq' %(args.prefix)
    cmds += 'rm ./nonhuman/%s_*human_bowtie2_contam.fastq; ' %(args.prefix)
    cmds += 'rm ./nonhuman/%s.log' %(args.prefix)

# step 6. remove vector sequences
# - input: ./nonhuman/%s.fastq (non-human reads)
# - output:
#   - ./nonhuman/%s.univec.txt: univec-contaminated sequence names
#   - ./nonhuman/%s.clean.fastq (univec-cleaned non-human reads)
if check_step(6, './nonhuman/%s.clean.fastq' %(args.prefix)):
    cmds += 'python /home/unix/csmillie/code/seq/fsq2fst.py ./nonhuman/%s.fastq | blastn -task blastn -reward 1 -penalty -5 -gapopen 3 -gapextend 3 -dust yes -soft_masking true -evalue 700 -searchsp 1750000000000 -outfmt 6 -db ~/aviv/db/univec/univec.fst | cut -f 1,7,8 > ./nonhuman/%s.univec.txt; ' %(args.prefix, args.prefix)
    cmds += 'python /home/unix/csmillie/code/seq/fasta_subset.py --fsq ./nonhuman/%s.fastq --remove %s.univec.txt --minlen 25 > ./nonhuman/%s.clean.fastq; ' %(args.prefix, args.prefix, args.prefix)


if args.p:
    print(cmds)
else:
    Submitter.submit([cmds])


# -------------------------------------------------------
# pipeline step 2: use kraken and blast to classify reads
# -------------------------------------------------------

if os.path.exists('./nonhuman/%s.clean.fastq' %(args.prefix)):
    cmd = 'rm ./nonhuman/%s.fastq' %(args.prefix)
    try:
        os.system(cmd)
    except:
        pass
else:
    quit()

cmds = []

# run kraken
if args.method in ['kraken', 'both']:

    cmd = ''
    
    # step 7. run kraken
    if check_step(7, '%s.kraken.out' %(args.prefix)):
        cmd += '/home/unix/csmillie/bin/kraken2/kraken2 --db %s ./nonhuman/%s.clean.fastq --report ./kraken/%s.kraken.report > ./kraken/%s.kraken.out; ' %(args.kraken_db, args.prefix, args.prefix, args.prefix)
    if check_step(7, '%s.bracken' %(args.prefix)):
        cmd += '/home/unix/csmillie/bin/bracken/bracken -d /home/unix/csmillie/aviv/db/kraken -i ./kraken/%s.kraken.report -o ./kraken/%s.bracken -r 100 -l G -t 5; ' %(args.prefix, args.prefix)
        cmds.append(cmd)

# run kneaddata
if args.method in ['knead', 'both']:
    orgs = 'archaea bacteria fungi protozoa viral'.split()
    for org in orgs:
        
        cmd = ''
        
        # step 8. run kneaddata and cleanup files
        # input: ./nonhuman/%s.clean.fastq (univec-cleaned non-human reads)
        # output: ./%s/%s_%s_bowtie2_contam.fastq (organism reads)
        if check_step(8, './%s/%s.%s.fastq' %(org, args.prefix, org)):
            cmd += 'kneaddata --input ./nonhuman/%s.clean.fastq -db /home/unix/csmillie/aviv/db/bowtie2/%s/%s --output %s --bypass-trim --output-prefix %s; ' %(args.prefix, org, org, org, args.prefix)
            cmd += 'rm ./%s/%s.fastq; ' %(org, args.prefix)
            cmd += 'mv ./%s/%s_%s_bowtie2_contam.fastq ./%s/%s.%s.fastq; ' %(org, args.prefix, org, org, args.prefix, org)
        
        # step 9. run blastn to map species
        # - require 95% identity, e-value <= 1e-4
        # input: ./%s/%s_%s_bowtie2_contam.fastq (organism reads)
        # output: ./%s/%s.%s.blastn (blastn output)
        if check_step(9, './%s/%s.%s.blastn' %(org, args.prefix, org)):
            cmd += "python /home/unix/csmillie/code/seq/fsq2fst.py ./%s/%s.%s.fastq | blastn -task blastn -db /home/unix/csmillie/aviv/db/refseq/v1/%s/%s.fna -outfmt 6 -evalue .0001 | awk '//{x[$1] += 1; if((x[$1] <= 10) && ($3 >= 95)){print $0}}' > ./%s/%s.%s.blastn; " %(org, args.prefix, org, org, org, org, args.prefix, org)
        
        # step 10. create blastn counts matrix
        if check_step(10, './%s/%s.%s.counts' %(org, args.prefix, org)):
            cmd += 'python /home/unix/csmillie/code/single_cell/microbe/2.blastn_counts.py --fastq ./%s/%s.%s.fastq --blast ./%s/%s.%s.blastn --sample %s --prefix %s --out ./%s/%s.%s.counts' %(org, args.prefix, org, org, args.prefix, org, args.prefix, org, org, args.prefix, org)
            cmds.append(cmd)

if args.p:
    print('\n'.join(cmds))
else:
    Submitter.submit(cmds)
