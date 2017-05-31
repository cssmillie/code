import argparse, re, subprocess

# read input arguments
parser = argparse.ArgumentParser()
parser.add_argument('--error', help='error file')
parser.add_argument('--regex', help='pattern')
args = parser.parse_args()

# find job and task id
jobs = {}
for line in open(args.error):
    if line.startswith('JobID'):
        if re.search(args.regex, line):
            line = line.rstrip().split('\t')
            job = line[1]
            task = line[2]
            if job not in jobs:
                jobs[job] = []
            jobs[job].append(task)

# print job information
print 'Found %s jobs and %s tasks' %(len(jobs), sum(map(len, jobs.values())))
for job in jobs:
    for task in jobs[job]:
        print 'Job %s, task %s' %(job, task)
        process = subprocess.Popen(['qacct -j %s -t %s' %(job, task)], stdout = subprocess.PIPE, shell=True)
        [out, error] = process.communicate()
        print out
        
