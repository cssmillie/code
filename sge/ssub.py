#!/usr/bin/env python

import argparse, os, os.path, re, stat, subprocess, sys, tempfile, time

'''
usage:
  cat list_of_commands.txt | ssub -n 100 -q short -m 8
  ssub -q long -m 8 "command1; command2; command3;"
'''

# set global variables
cluster = 'broad' # options = broad,amazon,coyote

# amazon header - sun grid engine
def amazon_header(n_jobs, error='error', queue='short', memory=None, array=False):
    h = '''#!/bin/bash
    source ~/.bashrc
    '''
    if array == True:
        h += '''
        #$ -t 1-%d
        #$ -j y
        #$ -o %s
        #$ -cwd
        ''' %(n_jobs, error, queue)
    h = re.sub('\n\s+', '\n', h)
    return h

# broad header - univa grid engine
def broad_header(n_jobs, error='error', queue='short', project=None, memory=None, array=False):
    h = '''#!/bin/bash
    source ~/.bashrc
    '''
    if array == True:
        h += '''
        #$ -t 1-%d
        #$ -j y
        #$ -o %s
        #$ -q %s
        #$ -cwd
        ''' %(n_jobs, error, queue)
        if project is not None:
            h += '\n#$ -P %s\n' %(project)
        if memory is not None:
            h += '\n#$ -l h_vmem=%dg\n' %(memory)
    h = re.sub('\n\s+', '\n', h)
    return h

# coyote header - torque
def coyote_header(n_jobs, error='error', queue='short', memory=None, array=False):
    h = '''#!/bin/bash
    source ~/.bashrc
    '''
    if array == True:
        h += '''
        #PBS -t 1-%d%%250
        #PBS -j oe
        #PBS -o %s
        #PBS -q %s
        cd $PBS_O_WORKDIR
        ''' %(n_jobs, error, queue)
        if memory is not None:
            h += '\n#PBS -l mem=%dgb' %(memory)
    h = re.sub('\n\s+', '\n', h)
    return h


class EmptyParser:
    def __init__(self):
        self.q = 'short'
        self.m = 0
        self.o = 'run'
        self.P = 'regevlab'
        self.H = ''
        self.p = False
        self.u = ''
        self.s = 'gold.broadinstitute.org'
        self.d = False
        self.commands = []


def parse_args():

    if __name__ == '__main__':
        
        # print usage statement
        usage = "\n\n  cat list_of_commands.txt | ssub -q long -m 8 -o test\n"
        usage +="  ssub -q short -m 16 'command 1; command 2; command 3;"
        
        # add command line arguments
        parser = argparse.ArgumentParser(usage = usage)
        parser.add_argument('-q', default='short', help='queue')
        parser.add_argument('-m', default=0, type=int, help='memory (gb)')
        parser.add_argument('-o', default='run', help='output prefix', required=True)
        parser.add_argument('-P', default='regevlab', help='project name')
        parser.add_argument('-H', default='', help='header lines')
        parser.add_argument('-p', default=False, action='store_true', help='print commands')
        parser.add_argument('-u', default='', help='username')
        parser.add_argument('-s', default='gold.broadinstitute.org', help='server')
        parser.add_argument('-d', default=False, action='store_true', help='direct submit (do not write scripts)')
        parser.add_argument('commands', nargs='?', default='')

        # parse arguments from stdin
        args = parser.parse_args()
        
    else:
        args = EmptyParser()
    
    return args


class Submitter():

    
    def __init__(self, cluster=cluster):
        
        # get command line arguments
        args = parse_args()
       
        # cluster parameters
        self.cluster = cluster
        self.m = args.m
        self.q = args.q
        self.o = args.o
        self.P = args.P
        self.p = args.p
        self.H = args.H
        self.u = args.u
        self.s = args.s
        self.d = args.d
        self.commands = args.commands
        
        if self.cluster == 'broad':
            self.header = broad_header
            self.submit_cmd = 'qsub'
            self.parse_job = lambda x: re.search('Your job-array (\d+)', x).group(1)
            self.stat_cmd = 'qstat -g d | egrep -v "^job|^-"'
            self.task_id = '$SGE_TASK_ID'
        
        if self.cluster == 'coyote':
            self.header = coyote_header
            self.submit_cmd = 'qsub'
            self.parse_job = lambda x: x.rstrip()
            self.stat_cmd = 'qstat -l'
            self.task_id = '$PBS_ARRAYID'

        if self.cluster == 'amazon':
            self.header = amazon_header
            self.submit_cmd = 'qsub'
            self.parse_job = lambda x: re.search('Your job-array (\d+)', x).group(1)
            self.stat_cmd = 'qstat -g d'
            self.task_id = '$SGE_TASK_ID'
        
        if self.m == 0:
            self.m = None
    

    def get_header(self, commands, error='error', array=False):
        h = self.header(n_jobs = len(commands), error=error, queue=self.q, project=self.P, memory=self.m, array=array)
        h = h + self.H + '\n'
        return h
    
    
    def system(self, commands):
        if self.u:
            cwd = os.path.realpath(os.getcwd())
            commands = ['ssh %s@%s "cd %s; %s"' %(self.u, self.s, cwd, command) for command in commands]    
        process = subprocess.Popen(commands, stdout=subprocess.PIPE, stderr=open(os.devnull, 'w'), shell=True)
        [out, err] = process.communicate()
        return out

    
    def qstat(self):
        # return job status
        commands = [self.stat_cmd]
        out = self.system(commands)
        out = [line for line in out.split('\n')]
        return out


    def njobs(self):
        out = self.qstat()
        njobs = len(out)
        return njobs
    
    
    def jobs_finished(self, job_ids):
        # check if jobs are finished
        out = self.qstat()
        for job in out:
            if job != '':
                job = job.strip().split()[0] # job is the first column of qstat
                for job_id in job_ids:
                    if job_id in job:
                        return False
        return True
    
    
    def qsub(self, fns):
        # submit jobs to the cluster
        job_ids = []
        for fn in fns:
            commands = ['%s %s' %(self.submit_cmd, fn)]
            out = self.system(commands)
            job_ids.append(self.parse_job(out))
            print('Submitting job %s' %(fn))
        return job_ids
    
    
    def write_array(self, commands):
        
        # get filenames
        fn1 = '.a.%s.sh' %(self.o)
        fn2 = '.a.%s.err' %(self.o)

        # warn if file exists
        if os.path.exists(fn1):
            print 'Warning: file %s exists' %(fn1)
        
        # write job array
        fh1 = open(fn1, 'w')
        fh1.write(self.get_header(commands=commands, error=fn2, array=True))
        for i, command in enumerate(commands):
            fh1.write('job_array[%d]="%s"\n' %(i+1, command))
        fh1.write('\neval ${job_array[%s]}\n\n' %(self.task_id))
        fh1.close()
        os.chmod(fn1, 0644)
        print('Writing array %s' %(fn1))
        
        # touch error log
        fh2 = open(fn2, 'a')
        fh2.close()
        
        # return array and error filenames
        return fn1, fn2
    
    
    def write_error(self, fn, commands, job_ids):
        # write header for output/error file
        error = open(fn, 'a')
        for i,command in enumerate(commands):
            error.write('\nJobID\t%s\t%s\t%s\n' %(job_ids[0], i+1, command))
        error.close()
    
    
    def submit(self, commands):
        # submit a job array to the cluster
        print '\nSubmitting jobs:'
        print '\t' + '\n\t'.join(commands)
        if len(commands) == 0:
            return []
        if type(commands) == str:
            commands = [commands]
        if self.p == False:
            if self.d == False:
                array_fn, error_fn = self.write_array(commands)
            else:
                if len(commands) > 1:
                    quit('Error: self.d == True and len(commands) > 1')
                array_fn = commands[0]
            job_ids = self.qsub([array_fn])
            self.write_error(error_fn, commands, job_ids)
            return job_ids
        elif self.p == True:
            print('\n'.join(commands))
            return []
    
    
    def wait(self, job_ids=[], max_jobs=''):
        # wait for job_ids to finish or number of jobs to go below max_jobs
        print '\nWaiting for jobs to finish:'
        print '\t' + '\n\t'.join(job_ids)
        if self.p == False:
            while True:
                time.sleep(5)
                if max_jobs:
                    if self.njobs() <= max_jobs:
                        break
                else:
                    if self.jobs_finished(job_ids):
                        break
    
    def submit_and_wait(self, commands, outs=[], retry=0, max_jobs=''):
        # submit job array and wait for it to finish
        queue = commands
        while retry >= 0:
            # submit jobs
            job_ids = self.submit(queue)
            self.wait(job_ids=job_ids, max_jobs=max_jobs)
            # check output
            queue = []
            if len(commands) == len(outs):
                for i in range(len(commands)):
                    if os.path.exists(outs[i]) and os.stat(outs[i]).st_size > 0:
                        continue
                    else:
                        queue.append(commands[i])
            if len(queue) == 0:
                return True
            retry -= 1
        return False

    
    def submit_pipeline(self, pipeline):
        # a pipeline is a list of lists of commands
        for commands in pipeline:
            self.submit_and_wait(commands)
    
    


def get_commands():
    ssub = Submitter(cluster)
    if ssub.commands != '':
        commands = [command.strip() for command in ssub.commands.split(';')]
    else:
        commands = [line.rstrip() for line in sys.stdin.readlines()]    
    return commands



if __name__ == '__main__':
    ssub = Submitter(cluster)
    commands = get_commands()
    ssub.submit(commands)

