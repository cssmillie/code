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
parser.add_argument('--fastq', help='FASTQ files', default='')
parser.add_argument('--bam', help='BAM files', default='')
parser.add_argument('--prefix', help='output prefix', default='')
parser.add_argument('--method', help='method', choices=['knead', 'kraken', 'both'], default='both')
parser.add_argument('--kraken_db', help='kraken database', default='/home/unix/csmillie/aviv/db/kraken')
parser.add_argument('--keep', help='keep all FASTQ files?', default=False, action='store_true')
parser.add_argument('-m', help='memory', type=int, default=64)
parser.add_argument('-t', help='runtime', default='8:00:00')
args = parser.parse_args()

# fix arguments
if args.fastq and args.bam:
    quit('error: must provide FASTQ or BAM')
if not args.fastq and not args.bam:
    quit('error: must provide FASTQ or BAM')

# get file lists
if args.fastq:
    args.fastq = args.fastq.split(',')
if args.bam:
    args.bam = args.bam.split(',')
if args.prefix:
    args.prefix = args.prefix.split(',')

# infer prefixes
if not args.prefix:
    if args.fastq:
        args.prefix = [os.path.basename(fi) for fi in args.fastq]
        args.prefix = [os.path.splitext(pi)[0] for pi in args.prefix]
    else:
        args.prefix = [os.path.basename(bi) for bi in args.bam]
        args.prefix = [os.path.splitext(pi)[0] for pi in args.prefix]


# --------------
# job management
# --------------

Submitter = ssub3.Submitter()
Submitter.m = [48,96,128]
Submitter.t = ['6:00:00','12:00:00','24:00:00']
Submitter.o = 'microbe'


# ------------------------------------
# pipeline step 1: get non-human reads
# ------------------------------------

for i in range(len(args.prefix)):

    # get prefix, bam, fastq
    # ----------------------
    pi = args.prefix[i]
    if args.bam:
        bi = args.bam[i]
    else:
        bi = None
    if args.fastq:
        fi = args.fastq[i]
    else:
        fi = None

    # ---------------------------------------------------------
    # part 1: get trimmed, complexity-filtered, non-human reads
    # ---------------------------------------------------------
    
    fn = './nonhuman/%s.clean.fastq' %(pi)
    if not os.path.exists(fn) or not os.path.getsize(fn) > 0:

        # 1. get unmapped reads from bam file
        # -----------------------------------
        # - input: bam file
        # - output: fastq file (*.fastq)
        if bi:
            fi = '%s.fastq' %(pi)
            cmd = 'python /home/unix/csmillie/code/seq/bam_to_fastq.py --bam %s --out %s --add_barcodes' %(bi, fi)
            Submitter.add_task(cmd, inputs=[bi], outputs=[fi])
        
        # 2. remove polyA tails and filter low complexity sequences
        # ---------------------------------------------------------
        # - input: %s.fastq (fastq file)
        # - output:
        #   - %s.good.tmp.fastq (trimmed and filtered sequences)
        #   - %s.bad.tmp.fastq (low complexity sequences)
        ii = fi
        print(pi)
        print(ii)
        oi = '%s.good.tmp.fastq' %(pi)
        cmd = 'prinseq-lite -fastq %s -lc_method dust -lc_threshold 16 -trim_tail_right 6 -trim_tail_left 6 -out_good %s.good.tmp -out_bad %s.bad.tmp' %(fi, pi, pi)
        Submitter.add_task(cmd, inputs=[ii], outputs=[oi])
        
        # 3. trim adapter sequences and require read length >= 20
        # -------------------------------------------------------
        # cleanup files to save space
        # - remove:
        #   - %s.fastq, %s.good.tmp.fastq, %s.bad.tmp.fastq
        # - keep
        #   - %s.trim.fastq (trimmed and length filtered sequences)
        ii = '%s.good.tmp.fastq' %(pi)
        oi = '%s.trim.fastq' %(pi)
        cmd = 'cutadapt -g "file:/home/unix/csmillie/bin/cutadapt.5p.TruSeq3-PE-2-SS2.fa" -m 20 -o %s.trim.fastq %s.good.tmp.fastq; ' %(pi, pi)
        if not args.keep:
            cmd += 'rm %s.fastq %s.good.tmp.fastq %s.bad.tmp.fastq;' %(pi, pi, pi)
        Submitter.add_task(cmd, inputs=[ii], outputs=[oi])
        
        # 4. remove human reads
        # ---------------------
        # - input: %s.trim.fastq (trimmed and filtered sequences)
        # - output:
        #   - ./nonhuman/%.fastq: non-human reads
        #   - ./nonhuman/%s_*human_bowtie2_contam.fastq: human reads
        #   - ./nonhuman/%s.log: logfile
        # cleanup files to save space:
        # - remove:
        #   - %s.trim.fastq, ./nonhuman/%s_*human_bowtie2_contam.fastq, ./nonhuman/%s.log
        # - keep:
        #   - ./nonhuman/%s.fastq: non-human reads
        ii = '%s.trim.fastq' %(pi)
        oi = './nonhuman/%s.fastq' %(pi)
        cmd = 'kneaddata --input %s.trim.fastq -db /home/unix/csmillie/aviv/db/bowtie2/human/human --bypass-trim --output nonhuman --output-prefix %s; ' %(pi, pi)
        if not args.keep:
            cmd += 'rm %s.trim.fastq ./nonhuman/%s_*human_bowtie2_contam.fastq ./nonhuman/%s.log; ' %(pi, pi, pi)
        Submitter.add_task(cmd, inputs=[ii], outputs=[oi])
        
        # 5. remove vector sequences
        # --------------------------
        # - input: ./nonhuman/%s.fastq (non-human reads)
        # - output:
        #   - ./nonhuman/%s.univec.txt: univec-contaminated sequence names
        #   - ./nonhuman/%s.clean.fastq (univec-cleaned non-human reads)
        ii = './nonhuman/%s.fastq' %(pi)
        oi = './nonhuman/%s.clean.fastq' %(pi)
        cmd = 'python /home/unix/csmillie/code/seq/fsq2fst.py ./nonhuman/%s.fastq | blastn -task blastn -reward 1 -penalty -5 -gapopen 3 -gapextend 3 -dust yes -soft_masking true -evalue 700 -searchsp 1750000000000 -outfmt 6 -db ~/aviv/db/univec/univec.fst | cut -f 1,7,8 > ./nonhuman/%s.univec.txt; ' %(pi, pi)
        cmd += 'python /home/unix/csmillie/code/seq/fasta_subset.py --fsq ./nonhuman/%s.fastq --remove ./nonhuman/%s.univec.txt --minlen 25 > ./nonhuman/%s.clean.fastq; ' %(pi, pi, pi)
        if not args.keep:
            cmd += 'rm ./nonhuman/%s.fastq; rm ./nonhuman/%s.univec.txt' %(pi, pi)
        Submitter.add_task(cmd, inputs=[ii], outputs=[oi])
        
    # ----------------------------------------------
    # part 2: use kraken and blast to classify reads
    # ----------------------------------------------
    
    # run kraken
    if args.method in ['kraken', 'both']:
        try:
            os.mkdir('kraken')
        except:
            pass
        
        # 6. run kraken
        ii = './nonhuman/%s.clean.fastq' %(pi)
        oi = './kraken/%s.kraken.out' %(pi)
        cmd = '/home/unix/csmillie/bin/kraken2/kraken2 --db %s ./nonhuman/%s.clean.fastq --report ./kraken/%s.kraken.report > ./kraken/%s.kraken.out; ' %(args.kraken_db, pi, pi, pi)
        Submitter.add_task(cmd, inputs=[ii], outputs=[oi])
            
        # 7. run bracken
        ii = './kraken/%s.kraken.report' %(pi)
        oi = './kraken/%s.bracken' %(pi)
        cmd = '/home/unix/csmillie/bin/bracken/bracken -d /home/unix/csmillie/aviv/db/kraken -i ./kraken/%s.kraken.report -o ./kraken/%s.bracken -r 100 -l G -t 5; ' %(pi, pi)
        Submitter.add_task(cmd, inputs=[ii], outputs=[oi])
        
    # run knead/blastn
    if args.method in ['knead', 'both']:
        orgs = 'archaea bacteria fungi protozoa viral'.split()
        for org in orgs:
            
            # 8. run kneaddata and cleanup files
            # ----------------------------------
            # input: ./nonhuman/%s.clean.fastq (univec-cleaned non-human reads)
            # output: ./%s/%s_%s_bowtie2_contam.fastq (organism reads)
            ii = './nonhuman/%s.clean.fastq' %(pi)
            oi = './%s/%s.%s.fastq' %(org, pi, org)
            cmd = 'kneaddata --input ./nonhuman/%s.clean.fastq -db /home/unix/csmillie/aviv/db/bowtie2/%s/%s --output %s --bypass-trim --output-prefix %s; ' %(pi, org, org, org, pi)
            cmd += 'mv ./%s/%s_%s_bowtie2_contam.fastq ./%s/%s.%s.fastq; ' %(org, pi, org, org, pi, org)
            if not args.keep:
                cmd += 'rm ./%s/%s.fastq; ' %(org, pi)            
            Submitter.add_task(cmd, inputs=[ii], outputs=[oi])
            
            # 9. run blastn to map species
            # ----------------------------
            # input: ./%s/%s_%s_bowtie2_contam.fastq (organism reads)
            # output: ./%s/%s.%s.blastn (blastn output)
            ii = './%s/%s.%s.fastq' %(org, pi, org)
            oi = './%s/%s.%s.blastn' %(org, pi, org)
            cmd = 'python /home/unix/csmillie/code/seq/fsq2fst.py ./%s/%s.%s.fastq | blastn -task blastn -db /home/unix/csmillie/aviv/db/refseq/v1/%s/%s.fna -outfmt 6 -evalue 10 > ./%s/%s.%s.blastn; ' %(org, pi, org, org, org, org, pi, org)
            Submitter.add_task(cmd, inputs=[ii], outputs=[oi])
            
            # 10. create blastn counts matrix
            # -------------------------------
            ii = './%s/%s.%s.blastn' %(org, pi, org)
            oi = ['./%s/%s.%s.sp_counts.txt' %(org, pi, org), './%s/%s.%s.gn_counts.txt' %(org, pi, org)]
            cmd = 'python /home/unix/csmillie/code/single_cell/microbe/2.blastn_counts.py --fastq ./%s/%s.%s.fastq --blast ./%s/%s.%s.blastn --sample %s --prefix %s --out ./%s/%s.%s' %(org, pi, org, org, pi, org, pi, org, org, pi, org)
            Submitter.add_task(cmd, inputs=[ii], outputs=oi)

Submitter.submit()
