import argparse, os, re, subprocess

me = 'csmillie'

def system(cmd, user='', p=False):
    # run system command and return output
    if user and user != me:
        cwd = os.path.realpath(os.getcwd())
        cmd = 'ssh %s@%s "cd %s; %s"' %(user, 'gold', cwd, cmd)
    if p == False:
        process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=open(os.devnull, 'w'), shell=True)
        [out, err] = process.communicate()
        return out
    else:
        print cmd


# ----------
# INPUT ARGS
# ----------

# command line arguments
parser = argparse.ArgumentParser()
parser.add_argument('-c', help='copies/duplicates', default=False, action='store_true')
parser.add_argument('-d', help='delete jobs (dry run)', default=False, action='store_true')
parser.add_argument('-D', help='delete jobs', default=False, action='store_true')
parser.add_argument('-e', help='error file', default=False, action='store_true')
parser.add_argument('-j', help='job id (default=all)', default='')
parser.add_argument('-q', help='queued', default=False, action='store_true')
parser.add_argument('-u', help='usernames', default='csmillie,mbiton,eugened,nrogel')
parser.add_argument('-x', help='job regex', default='')
args = parser.parse_args()

# process arguments
users = args.u.split(',')


# --------------------
# RUNNING JOBS (QSTAT)
# --------------------

# get qstat output
cmd = 'qstat -g d -u %s | grep -v QRLOGIN | egrep -v "^job|^-"' %(','.join(users))
out = system(cmd).split('\n')

# map jobs to tasks
x = {}
for line in out:
    if args.q and 'qw' not in line:
        continue
    line = line.rstrip().split()
    if len(line) == 0:
        continue
    job = line[0]
    if args.j and args.j != job:
        continue
    task = line[-1]
    x[job] = x.get(job, []) + [task]


# ----------------------
# JOB DETAILS (QSTAT -J)
# ----------------------


def get_time(line):
    h,m,s = map(float, re.search('wallclock=(\d+):(\d+):(\d+),', line).groups())
    return 60*h + m + s/60.


# aggregate jobs by commands
y = {}
for job in x:
    
    try:
        
        cmd = 'qstat -j %s' %(job)
        out = system(cmd).split('\n')
        
        # parse details
        cwd = [line.rstrip().split()[-1] for line in out if line.startswith('cwd')][0]
        arr = [line.rstrip().split()[-1] for line in out if line.startswith('script_file:')][0]
        err = [line.rstrip().split()[-1] for line in out if line.startswith('stdout_path_list')][0]
        mem = [line.rstrip().split()[-1] for line in out if line.startswith('hard resource_list:')][0]
        user = [line.rstrip().split()[-1] for line in out if line.startswith('owner:')][0]
        time = [get_time(line) for line in out if 'wallclock' in line]
        
        # fix vars
        arr = os.path.join(cwd, arr)
        err = re.sub('NONE:', '', err)
        err = os.path.join(cwd, err)
        if len(time) == 0:
            time = 0
        else:
            time = time[0]
        
        # job commands
        for line in open(arr):
            m = re.search('array\[(\d+)\]="(.*)"', line)
            if m:
                task = m.group(1)
                command = m.group(2)
                if task in x[job]:
                    if args.x and not re.search(args.x, command):
                        continue
                    if command not in y:
                        y[command] = []
                    y[command].append([user, job, task, mem, time, command, err])
    except:
        print 'Could not map %s' %(job)


# ----------
# DUPLICATES
# ----------

if args.c:
    for command in y:
        max_time = max([yi[4] for yi in y[command]])
        y[command] = [yi for yi in y[command] if yi[4] < max_time]


# ------
# OUTPUT
# ------    

for command in y:
    for [user, job, task, mem, time, command, err] in y[command]:
        if args.d:
            system('qdel %s -t %s' %(job, task), user=user, p=True)
        elif args.D:
            system('qdel %s -t %s' %(job, task), user=user, p=False)
        elif args.e:
            print err
        else:
            print '\t'.join(map(str, [user, job, task, mem, time, command]))
