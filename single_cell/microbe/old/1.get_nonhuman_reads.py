import argparse, re, ssub

bams = []
sids = []

for line in open('/home/unix/csmillie/Gut_Human/csmillie/data/Sample_Info.txt'):
    line = line.rstrip().split()
    folder = line[0]
    subject = line[1]
    if subject == 'N5':
        continue
    for site in 'Epi_A Epi_B LP_A LP_B'.split():
        bams.append('/home/unix/csmillie/Gut_Human/data/%s/%s/outs/possorted_genome_bam.bam' %(folder, site))
        sids.append('%s.%s' %(subject, re.sub('_', '', site)))

# Setup ssub pipeline
Submitter = ssub.Submitter()
Submitter.m = 64
Submitter.q = 'long'
Submitter.o = 'nonhuman'
Submitter.p = False
pipeline = []

# Get unmapped reads, trim polyA tails, and remove low complexity reads
cmds = []
for i in range(len(bams)):
    
    bam = bams[i]
    out = sids[i]
    
    # The bam2fq command creates the following files:
    # ./%s.good.tmp.fastq
    # ./%s.bad.tmp.fastq
    # Remove the bad sequences to save space
    # rm %s.bad.tmp.fastq
    # Move the good sequences to simplify names
    # mv %s.good.tmp.fastq %s.fastq
    cmdi = 'samtools bam2fq -f 4 %s > %s.fastq; prinseq-lite -fastq %s.fastq -lc_method dust -lc_threshold 16 -trim_tail_right 6 -trim_tail_left 6 -out_good %s.good.tmp -out_bad %s.bad.tmp; rm %s.bad.tmp.fastq; mv %s.good.tmp.fastq %s.fastq' %(bam, out, out, out, out, out, out, out)
    
    # Remove human reads
    cmdj = 'kneaddata --input %s.fastq -db /home/unix/csmillie/aviv/db/bowtie2/human/human --output human --output-prefix %s' %(out, out)
    
    # This command creates the following output files
    # ./human/%s.fastq
    # ./human/%s.log
    # ./human/%s.trimmed.fastq
    # ./human/%s_human_bowtie2_contam.fastq
    
    # Remove the trimmed FASTQ file to save space
    cmdk = 'rm ./human/%s.trimmed.fastq' %(out)

    cmds.append('%s; %s; %s' %(cmdi, cmdj, cmdk))

pipeline.append(cmds)


# Find non-human reads
cmds = []
orgs = 'archaea bacteria fungi protozoa viral'.split()

for i in range(len(bams)):
    
    bam = bams[i]
    out = sids[i]
    
    for org in orgs:
        cmdi = 'kneaddata --input ./human/%s.fastq -db /home/unix/csmillie/aviv/db/bowtie2/%s/%s --output %s --output-prefix %s' %(out, org, org, org, out)
        cmdj = 'rm ./%s/%s.fastq' %(org, out)
        cmdk = 'rm ./%s/%s.trimmed.fastq' %(org, out)
        cmds.append('%s; %s; %s' %(cmdi, cmdj, cmdk))

pipeline.append(cmds)

# Submit pipeline
Submitter.submit_pipeline(pipeline)
