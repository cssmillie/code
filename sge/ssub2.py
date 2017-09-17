#!/usr/bin/env python
import argparse, datetime, itertools, os, os.path, random, re, stat, subprocess, sys, tempfile, time


# ---------------
# input arguments
# ---------------


def parse_args():

    # add command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('-m', default=0, type=int, help='minimum memory (gb)')
    parser.add_argument('-M', default=0, type=int, help='maximum memory (gb)')        
    parser.add_argument('-q', default='short', help='queue')
    parser.add_argument('-Q', default='long', help='retry queue')
    parser.add_argument('-o', default='run', help='output prefix', required=True)
    parser.add_argument('-P', default='regevlab', help='project name')
    parser.add_argument('-H', default='', help='header lines')
    parser.add_argument('-u', default='', help='username')
    parser.add_argument('-s', default='gold.broadinstitute.org', help='server')
    parser.add_argument('-p', default=False, action='store_true', help='print commands')        
    parser.add_argument('-d', default=False, action='store_true', help='direct submit')
    parser.add_argument('-w', default=600, type=float, help='pipeline wait time (sec)')
    parser.add_argument('-r', default=0, type=int, help='max retry')
    parser.add_argument('-W', default=float('inf'), type=float, help='max inactivity (sec)')
    parser.add_argument('-f', default='.', help='folder for sh and err files')
    parser.add_argument('commands', nargs='?', default='')
    
    if __name__ == '__main__':
        # parse arguments
        args = parser.parse_args()
    else:
        # use defaults
        args = parser.parse_args(['-o', 'test'])
    
    return args



# ---------------
# main subroutine
# ---------------


if __name__ == '__main__':

    # initialize new task submitter
    submitter = Submitter()

    # add commands from stdin
    for line in sys.stdin.readlines():
        command = line.rstrip()
        submitter.add_task(command)

    # submit jobs
    submitter.submit()



# ------------------
# global definitions
# ------------------


def check_files(fns):
    for fn in fns:
        if not os.path.exists(fn) or os.stat(fn).st_size == 0:
            return False
    return True


def get_filename(prefix, suffix='txt', path='.'):
    # make directory
    if not os.path.exists(path):
        os.mkdir(path)
    # find filenames
    i = 0
    while True:
        i += 1
        fn = '%s/%s.%s.%s' %(path, prefix, i, suffix)
        if not os.path.exists(fn):
            break
    return fn



# ---------
# write log
# ---------


class Log():

    def __init__(self, prefix='.a.run', suffix='log', path='.'):
        self.start = datetime.datetime.now()
        self.out = get_filename(prefix=prefix, suffix=suffix, path=path)
        self.out = open(self.out, 'w')
    
    def write(self, text):
        text = re.sub('^(?=[^\s])', '    ', text)
        text = re.sub('\n', '\n    ', text)
        time = datetime.datetime.now()
        sec = (time - self.start).total_seconds()
        text = '[%s | %.2f min | %d sec]\n%s\n' %((time).strftime('%Y-%m-%d %H:%M:%S'), sec/60., sec, text)
        print text,
        self.out.write(text)



# -----------
# task object
# -----------


class Task():
    
    def __init__(self, command, infiles=[], outfiles=[], intasks=[], m=[0], q=['short'], u=['']):
        
        self.command = command

        # each task requests separate resources
        self.m = None # memory
        self.q = None # queue
        self.u = None # username
        
        # lists of resources to iterate through
        fix_type = lambda x: x if type(x) == list else [x]
        self.M = fix_type(m)[:]
        self.Q = fix_type(q)[:]
        self.U = fix_type(u)[:]
        
        # job scheduling parameters
        self.infiles = infiles
        self.intasks = intasks
        self.outfiles = outfiles
        self.status = 'waiting'
        self.uid = ''
        self.n = 0
        
        # set initial resources to M[0], Q[0], U[0]
        self = self.update_resources()
    
    
    def update_resources(self, job_id='', task_id=''):
        
        # get next parameters in M, Q, and U
        if len(self.M) > 0:
            self.m = self.M.pop(0)
        if len(self.Q) > 0:
            self.q = self.Q.pop(0)
        if len(self.U) > 0:
            self.u = self.U.pop(0)
        
        # update uid and increment n
        if job_id != '' and task_id != '':
            self.uid = '%s;%s' %(job_id, task_id)
            self.n += 1
        
        # return task with increased resources
        return self
    
    
    def __repr__(self):
        x = vars(self)
        return '\n' + self.__class__.__name__ + '\n' + '\n'.join([k.ljust(10) + str(x[k]) for k in sorted(x)]) + '\n'

    
    def check_infiles(self):
        # test if all infiles exist
        return check_files(self.infiles)

    
    def check_outfiles(self):
        # test if all outfiles exist
        return check_files(self.outfiles)

    
    def check_intasks(self):
        # test if all intasks finished
        for task in self.intasks:
            if task.status != 'finished':
                return False
        return True
    
    
    def job_id(self):
        # return job id
        return self.uid.split(';')[0]

    
    def task_id(self):
        # return task id
        return self.uid.split(';')[1]

    
    def resources(self):
        # return resources string
        return '%s.%s.%s' %(self.u, self.q, self.m)

    
# --------------
# task submitter
# --------------


class Submitter():
    
    
    def __init__(self):
        
        # command line arguments
        args = parse_args()
        
        # submission parameters
        self.tasks = [] # tasks
        self.m = args.m # memory
        self.q = args.q # queue
        self.o = args.o # output
        self.P = args.P # project
        self.H = args.H # header
        self.u = args.u # user
        self.s = args.s # server
        self.p = args.p # print
        self.d = args.d # direct
        self.w = args.w # wait time
        self.r = args.r # max retry        
        self.W = args.W # max inactivity
        self.f = args.f # folder
        self.inactivity = time.time()
        
        # fix arguments
        self.m = map(int, str(self.m).split(','))
        self.q = self.q.split(',')
        self.u = self.u.split(',')
        
        # account setup
        self.me = 'csmillie'
        self.users = ['', 'eugened', 'mbiton']
        
        # cluster parameters
        self.stat_cmd = 'qstat -g d -u %s | egrep -v "^job|^-"' %(','.join(self.users + [self.me]))
        self.submit_cmd = 'qsub'
        self.parse_job = lambda x: re.search('Your job-array (\d+)', x).group(1)
        
        # array writer
        self.log = Log()
        self.writer = Writer()
        
        # initialize job objects
        for command in args.commands:
            self.add_task(command)
    
    
    def write_log(self, text):
        self.log.write(text)

    
    def system(self, cmd, user=''):
        # run system command and return output
        if user:
            cwd = os.path.realpath(os.getcwd())
            cmd = 'ssh %s@%s "cd %s; %s"' %(user, self.s, cwd, cmd)
        process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=open(os.devnull, 'w'), shell=True)
        [out, err] = process.communicate()
        return out
    
    
    def qstat(self):
        # job status
        out = self.system([self.stat_cmd])
        out = [line for line in out.split('\n')]
        return out    

    
    def qsub(self, fn, user='', task_ids=[]):
        # submit jobs
        if len(task_ids) == 0:
            out = self.system('%s %s' %(self.submit_cmd, fn), user=user)
        else:
            out = self.system('%s -t %s %s' %(self.submit_cmd, ','.join(task_ids), fn), user=user)
        job_id = self.parse_job(out)
        return job_id
        
    
    def new_task(self, command='', infiles=[], outfiles=[], intasks=[], m=None, q=None, u=None):
        assert type(command) == str
        if m is None:
            m = self.m
        if q is None:
            q = self.q
        if u is None:
            u = self.u
        task = Task(command, infiles=infiles, outfiles=outfiles, intasks=intasks, m=m, q=q, u=u)
        return task
    
    
    def add_task(self, task):
        assert task.__class__.__name__ == 'Task'
        self.write_log('Task "%s" (m=%s, q=%s, u=%s)' %(task.command, task.m, task.q, task.u))
        self.tasks.append(task)
    
    
    def add_command(self, command, infiles=[], outfiles=[], intasks=[], m=None, q=None, u=None):
        assert type(command) == str
        task = self.new_task(command, infiles=infiles, outfiles=outfiles, intasks=intasks, m=m, q=q, u=u)
        self.add_task(task)


    def get_running_uids(self):
        x = {}
        for line in self.qstat():
            line = line.rstrip().split()
            if len(line) > 0:
                uid = '%s;%s' %(line[0], line[-1])
                x[uid] = 1
        return x
    
    
    def update_tasks(self):
        
        # single tasks
        for task in self.tasks:
            if not task.check_infiles():
                task.status = 'waiting'
            elif not task.uid:
                task.status = 'ready'
            elif task.uid in self.get_running_uids():
                task.status = 'running'
            elif not task.check_outfiles():
                if task.n <= self.r:
                    task.status = 'retry'
                else:
                    task.status = 'failed'
            else:
                task.status = 'finished'        
        
        # pipelines
        for task in self.tasks:
            if not task.check_intasks():
                task.status = 'waiting'
    
    
    def group_tasks_by_status(self, tasks):
        # split tasks into status groups
        self.update_tasks()
        x = {}
        for status in ['waiting', 'ready', 'running', 'retry', 'failed', 'finished']:
            x[status] = []
        for task in self.tasks:
            x[task.status] = x.get(task.status, []) + [task]
        return x
    
    
    def group_tasks_by_resources(self, tasks):
        # split tasks into resource groups
        self.update_tasks()
        x = {}
        for task in tasks:
            resource = task.resources()
            x[resource] = x.get(resource, []) + [task]
        return x
    
    
    def assign_tasks_to_users(self, tasks, users, shuffle=True):
        random.shuffle(users)
        users = itertools.cycle(users)
        if shuffle == True:
            random.shuffle(tasks)
        for task in tasks:
            task.u = next(users)
    
    
    def summarize_tasks(self):

        # summarize status
        x = self.group_tasks_by_status(self.tasks)
        k = 'waiting ready running retry failed finished'.split()
        text =  ''.join([ki.ljust(10) for ki in k]) + '\n'
        text += ''.join([str(len(x[ki])).ljust(10) for ki in k]) + '\n'

        # summarize resources
        x = self.group_tasks_by_resources(self.tasks)
        text += ''.join([ki.ljust(20) for ki in x]) + '\n'
        text += ''.join([str(len(x[ki])).ljust(20) for ki in x])
        self.write_log(text)
    
    
    def submit(self, commands=[], tasks=[], status='ready', randomize=True):
        
        # add commands
        for command in commands:
            self.add_command(command)

        # add tasks
        for task in tasks:
            self.add_task(task)
        
        # assign tasks to users
        if randomize:
            self.assign_tasks_to_users(self.tasks, self.users[:], shuffle=True)
        
        # get tasks
        tasks = self.group_tasks_by_status(self.tasks)[status]
        
        # split into resource groups
        tasks = self.group_tasks_by_resources(tasks)
        
        # submit each resource group
        for group in tasks:
            
            # get tasks and commands
            ti = tasks[group]
            ci = [task.command for task in ti if task.n <= self.r]
            mi = ti[0].m # memory
            qi = ti[0].q # queue
            ui = ti[0].u # username
            
            if len(ci) == 0:
                return []
            
            # case 1: print submit
            if self.p == True:
                return []
            
            # case 2: direct submit
            elif self.d == True:
                array_fn = ci[0]
                job_id = self.qsub(array_fn, user=self.u)
            
            # case 3: normal submit
            else:
                array_fn, error_fn = self.writer.write_array(ci, prefix=self.o, q=qi, P=self.P, m=mi, H=self.H, path=self.f)
                job_id = self.qsub(array_fn, user=ui)
            
            self.write_log('submitting array %s\n(user: %s, queue: %s, mem: %s)\n%s' %(array_fn, ui, qi, mi, '\n'.join(ci)))
            self.writer.write_error(error_fn, commands, job_id)
            
            # update objects
            for i,tj in enumerate(ti):
                tj.update_resources(job_id=job_id, task_id=i+1)
    
    
    def check_finished(self):
        tasks = self.group_tasks_by_status(self.tasks)
        
        # case 1: all tasks finished
        if len(tasks['finished']) + len(tasks['failed']) == len(self.tasks):
            self.write_log('all tasks finished')
            return True
        
        # case 2: inactivity time out
        if len(tasks['running']) > 0:
            self.inactivity = time.time()
        ti = time.time() - self.inactivity
        if ti > self.W:
            self.write_log('time out: %.2f sec' %(ti))
            return True
        
        # case 3: tasks not finished
        return False
    
    
    def submit_pipeline(self, commands=[], tasks=[], randomize=True):
        
        # add tasks to submitter
        for command in commands:
            self.add_task(command)
        
        # run jobs until finished
        while not self.check_finished():
            # summarize tasks
            self.summarize_tasks()
            # submit 'ready' jobs
            self.submit(status='ready', randomize=randomize)
            # submit 'retry' jobs (increased resources)
            self.submit(status='retry', randomize=randomize)
            # wait until next submission
            time.sleep(self.w)
    


# -----------------
# write task arrays
# -----------------


class Writer():
    
    
    def __init__(self):
        self.task = '$SGE_TASK_ID'
    
    
    def get_header(self, commands, error='error', q='short', P=None, m=0, H='', array=False):
        
        h = '''
        #!/bin/bash
        source ~/.bashrc
        '''
        
        if array == True:
            h += '''
            #$ -t 1-%d
            #$ -j y
            #$ -o %s
            #$ -q %s
            #$ -cwd
            ''' %(len(commands), error, q)
        
        if P is not None:
            h += '''
            #$ -P %s
            ''' %(P)
        
        if m != 0:
            h += '''
            #$ -l h_vmem=%dg
            ''' %(m)
        
        h = re.sub('\n\s+', '\n', h)
        h = '%s%s\n' %(h, H)
        return h
        
    
    def write_array(self, commands, prefix, q='short', P=None, m=0, H='', path='.'):
        
        # get filenames
        prefix = '.a.%s' %(prefix)
        fn1 = get_filename(prefix=prefix, suffix='sh', path=path)
        fn2 = re.sub('.sh$', '.err', fn1)
        
        # write job array
        fh = open(fn1, 'w')
        fh.write(self.get_header(commands=commands, q=q, P=P, m=m, H=H, error=fn2, array=True))
        for i,command in enumerate(commands):
            fh.write('job_array[%d]="%s"\n' %(i+1, command))
        fh.write('\neval ${job_array[%s]}\n\n' %(self.task))
        fh.close()
        os.chmod(fn1, 0644)
        
        # touch error log
        fh = open(fn2, 'a')
        fh.close()
        
        # return array and error filenames
        return fn1, fn2
    
    
    def write_error(self, error, commands, job_id):
        # write header for error file
        error = open(error, 'a')
        for i,command in enumerate(commands):
            error.write('\nJobID\t%s\t%s\t%s\n' %(job_id, i+1, command))
        error.close()
