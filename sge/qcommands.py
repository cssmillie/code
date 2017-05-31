import argparse, re, subprocess

uname = 'csmillie'

# Get command line arguments
parser = argparse.ArgumentParser()
parser.add_argument('-j', help='job id (default=all)', default='')
args = parser.parse_args()

# Get qstat output
cmd = 'qstat'
out = subprocess.check_output(cmd).split('\n')

# Map jobs to tasks
job2task = {}
for line in out:
    if re.search(uname, line):
        try:
            line = line.rstrip().split()
            job_id = line[0]
            task_id = line[-1]
            if args.j and job_id != args.j:
                continue
            if job_id not in job2task:
                job2task[job_id] = []
            if ':' not in task_id:
                job2task[job_id].append(task_id)
            else:
                [beg,end] = re.search('(\d+)-(\d+)', task_id).groups()
                for i in range(int(beg), int(end)+1):
                    job2task[job_id].append(str(i))
        except:
            continue

# Get commmands
for job_id in job2task:
    cmd = ['qstat', '-j', job_id]
    out = subprocess.check_output(cmd).split('\n')
    afn = [line.strip().split()[-1] for line in out if 'submit_cmd' in line][0]
    mem = [line.strip().split()[-1] for line in out if 'h_vmem' in line][0]
    jfn = re.sub('\.a\.', '.j.', afn)
    try:
        for line in open(jfn):
            m = re.search('job_array\[(\d+)\]="(.*)"', line)
            if m:
                task_id = m.group(1)
                if task_id in job2task[job_id]:
                    print job_id, mem, m.group(2)
    except:
        continue
