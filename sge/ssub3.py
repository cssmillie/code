#!/usr/bin/env python
import argparse, datetime, itertools, os, os.path, random, re, stat, subprocess, sys, tempfile, time


# ---------------
# input arguments
# ---------------


def parse_args():

    # add command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('-m', default=8, type=int, help='memory (gb)')
    parser.add_argument('-q', default='short', help='queue')
    parser.add_argument('-o', default='run', help='output prefix', required=True)
    parser.add_argument('-P', default='regevlab', help='project name')
    parser.add_argument('-H', default='', help='header lines')
    parser.add_argument('-u', default='', help='username')
    parser.add_argument('-s', default='gold.broadinstitute.org', help='server')
    parser.add_argument('-p', default=False, action='store_true', help='print commands')        
    parser.add_argument('-d', default=False, action='store_true', help='direct submit')
    parser.add_argument('-w', default=600, type=float, help='pipeline wait time (sec)')
    parser.add_argument('-r', default=0, type=int, help='max retry')
    parser.add_argument('-R', default=True, action='store_false', help='do not randomize users')
    parser.add_argument('-W', default=float('inf'), type=float, help='max inactivity (sec)')
    parser.add_argument('-g', default=False, action='store_true', help='group tasks')
    parser.add_argument('-v', default=True, action='store_false', help='quiet')
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


# -----------
# task object
# -----------


class Task():

    
    def __init__(self, command, inputs=[], outputs=[], m=[0], q=['short'], u=['']):
        
        # job tracking
        self.command = command
        self.array = None
        self.error = None
        self.job_id = None
        self.task_id = None
        self.n = -1
        
        # job scheduling
        self.inputs = inputs
        self.outputs = outputs
        self.status = 'waiting'
        
        # current resources
        self.m = None # memory
        self.q = None # queue
        self.u = None # username
        
        # list of resources
        fix_type = lambda x: x if type(x) == list else [x]
        self.M = fix_type(m)[:]
        self.Q = fix_type(q)[:]
        self.U = fix_type(u)[:]

        # initialize resources
        self = self.add_iter(job_id=None)
    
    
    def __repr__(self):
        x = vars(self)
        return '\n' + self.__class__.__name__ + '\n' + '\n'.join(['%s:'.ljust(10) %(k) + str(x[k]) for k in sorted(x)]) + '\n'
    
    
    def add_iter(self, job_id=None):
        self.job_id = job_id
        self.n += 1
        self.m = self.m if len(self.M) == 0 else self.M.pop(0)
        self.q = self.q if len(self.Q) == 0 else self.Q.pop(0)
        self.u = self.u if len(self.U) == 0 else self.U.pop(0)
        return self
    
    
    def check_inputs(self):
        return check_files(self.inputs)
    
    
    def check_outputs(self):
        if len(self.outputs) == 0:
            return False
        return check_files(self.outputs)
    
    
    def uid(self):
        return '%s;%s' %(self.job_id, self.task_id)
    
    
    def resources(self):
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
        self.w = args.w # wait time
        self.r = args.r # max retry
        self.R = args.R # randomize users
        self.W = args.W # max inactivity
        self.g = args.g # group tasks
        self.v = args.v # verbose
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
        
        # initialize job objects
        for command in args.commands:
            self.add_task(command)
    
    
    def write_log(self, text):
        self.writer.write_log(text, verbose=self.v)
    
    
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
    
    
    def qsub(self, tasks):
        
        # qsub parameters
        commands = [task.command for task in tasks]
        def f(x):
            x = list(set([getattr(task, x) for task in tasks]))
            assert len(x) == 1
            return x[0]
        array = f('array')
        error = f('error')
        m = f('m')
        q = f('q')
        u = f('u')
        t = ','.join(map(str, sorted([task.task_id for task in tasks])))        
        
        # submit job
        cmd = '%s -o %s -j y -l h_vmem=%sg -q %s -t %s -P %s %s' %(self.submit_cmd, error, m, q, t, self.P, array)
        out = self.system(cmd, user=u)
        
        # update tasks
        job_id = self.parse_job(out)
        for task in tasks:
            task.add_iter(job_id=job_id)
        
        # write log
        self.write_log('Job %s: "%s" (u=%s)\n%s' %(job_id, cmd, u, '\n'.join(commands)))
        
    
    def add_task(self, command='', inputs=[], outputs=[], m=None, q=None, u=None):
        assert type(command) == str
        
        # use default resources
        m = self.m if m is None else m
        q = self.q if q is None else q
        u = self.u if u is None else u
        
        # create and return task
        task = Task(command, inputs=inputs, outputs=outputs, m=m, q=q, u=u)
        self.tasks.append(task)
    
    
    def running_uids(self):
        x = {}
        for line in self.qstat():
            line = line.rstrip().split()
            if len(line) > 0:
                uid = '%s;%s' %(line[0], line[-1])
                x[uid] = 1
        return x
    
    
    def update_tasks(self):
        # update task status with single qstat
        for task in self.tasks:
            if task.uid() in self.running_uids():
                task.status = 'running'
            elif task.check_outputs():
                task.status = 'finished'
            elif not task.check_inputs():
                task.status = 'waiting'
            elif task.n <= self.r:
                task.status = 'ready'
            else:
                task.status = 'failed'
        return self
    
    
    def group_tasks_by_status(self, tasks):
        # split tasks into status groups
        self = self.update_tasks()
        x = {}
        for status in ['waiting', 'ready', 'running', 'failed', 'finished']:
            x[status] = []
        for task in self.tasks:
            x[task.status] = x.get(task.status, []) + [task]
        return x
    
    
    def group_tasks_by_resources(self, tasks):
        # split tasks into resource groups
        self = self.update_tasks()
        x = {}
        for task in tasks:
            resources = task.resources()
            x[resources] = x.get(resources, []) + [task]
        return x
    
    
    def summarize_tasks(self):
        
        # summarize status
        order = ['waiting', 'ready', 'running', 'failed', 'finished']
        x = self.group_tasks_by_status(self.tasks)
        u = ['%s: %s' %(ki, len(x[ki])) for ki in order]
        
        # summarize resources
        x = self.group_tasks_by_resources(self.tasks)
        v = ['%s: %s' %(ki, len(x[ki])) for ki in x]
        
        # write to log
        text = ', '.join(u) + '\n' + ', '.join(v)
        self.write_log(text)
    
    
    def submit(self, x=[], wait=True, out=None, group=None):
        
        # add tasks/commands
        for xi in x:
            self.add_task(xi)
        self.o = self.o if out is None else out
        self.g = self.g if group is None else group
        
        if len(self.tasks) == 0:
            return []
        
        # initialize log
        self.writer = Writer(prefix=self.o, path=self.o, verbose=self.v)
        
        # write array
        array = self.writer.write_array(tasks=self.tasks, H=self.H)

        # randomize users
        if self.R:
            for task in self.tasks:
                task.u = random.choice(self.users)
        
        while not self.check_finished():
            
            # summarize all tasks
            self.summarize_tasks()            
            
            # get 'ready' tasks
            tasks = self.group_tasks_by_status(self.tasks)['ready']

            # split into resource groups
            if self.g == True:
                groups = self.group_tasks_by_resources(tasks)
            else:
                groups = [[task] for task in tasks]
            
            # submit each resource group
            for tasks in groups:
                self.writer.write_error(tasks=tasks)
                self.qsub(tasks)
            
            # go to sleep
            if wait == False:
                break
            self.writer.log.flush()
            time.sleep(self.w)
        
        self.writer.log.close()
    
    
    def check_finished(self):

        # get task groups
        tasks = self.group_tasks_by_status(self.tasks)

        # case 1: finished
        if len(tasks['finished']) + len(tasks['failed']) == len(self.tasks):
            self.write_log('finished')
            return True
        
        # case 2: inactivity
        if len(tasks['running']) > 0:
            self.inactivity = time.time()
        ti = time.time() - self.inactivity
        if ti > self.W:
            self.write_log('time out: %.2f sec' %(ti))
            return True
        
        # case 3: not finished
        return False
    

# -----------
# write files
# -----------


class Writer():
    
    
    def __init__(self, prefix, path, verbose=True):
        
        # filename info
        self.prefix = prefix
        self.path = path

        # initialize log
        self.start = datetime.datetime.now()
        self.log = open(get_filename(prefix=prefix, suffix='log', path=path), 'w')
        self.verbose = verbose

    
    def array_header(self, H=''):
        
        # write header
        h = '''
        #!/bin/bash
        source ~/.bashrc
        #$ -cwd
        \n''' + H
        
        # fix whitespace
        h = re.sub('\n\s+', '\n', h)
        
        return h
    
    
    def write_array(self, tasks, H=''):
        
        # get filename
        array = get_filename(prefix=self.prefix, suffix='sh', path=self.path)
        
        # array text
        text = self.array_header(H) + '\n'
        for i,task in enumerate(tasks):
            text += 'array[%d]="%s"\n' %(i+1, task.command)
        text += '\neval ${array[$SGE_TASK_ID]}\n\n'
        
        # write array
        with open(array, 'w') as fh:
            fh.write(text)
        os.chmod(array, 0644)
        
        # update tasks
        for i,task in enumerate(tasks):
            task.array = array
            task.task_id = i+1
        
        # return filename
        return array
    
    
    def write_error(self, tasks):
        
        # get prefix
        array = list(set([task.array for task in tasks]))
        assert len(array) == 1
        array = array[0]
        prefix = re.sub('.*\/', '', re.sub('.sh$', '', array))
        
        # get filename
        if len(tasks) > 1:
            error = get_filename(prefix='.'+prefix, suffix='err', path=self.path)
        else:
            error = re.sub(r'^(.*\/)', r'\1.', re.sub('.sh$', '.%s.err' %(tasks[0].task_id), array))
        
        # error text
        text = ''
        for task in tasks:
            text = 'Task %d: "%s" (m=%s, q=%s, u=%s)' %(task.task_id, task.command, task.m, task.q, task.m) + '\n'
        
        # write error
        with open(error, 'a') as fh:
            fh.write(text)
        
        # update tasks
        for task in tasks:
            task.error = error
        
        # return filename
        return error
    
    
    def write_log(self, text, verbose=True):
        
        # log text
        text = re.sub('^(?=[^\s])', '    ', text)
        text = re.sub('\n', '\n    ', text)
        time = datetime.datetime.now()
        sec = (time - self.start).total_seconds()
        text = '[%s | %.2f min | %d sec]\n%s\n' %((time).strftime('%Y-%m-%d %H:%M:%S'), sec/60., sec, text)
        
        # write log
        if self.verbose:
            print text,
        self.log.write(text)
