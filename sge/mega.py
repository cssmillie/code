import argparse, os, re, subprocess


def system(cmd, user='csmillie', p=False):
    # run system command and return output
    process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=open(os.devnull, 'w'), shell=True)
    [out, err] = process.communicate()
    return out

# get running jobs
cmd = 'qstat -g d | grep -v QRLOGIN | egrep -v "^job|^-"'
out = system(cmd).split('\n')
out = [line.rstrip().split() for line in out]
jobs = [line[0] for line in out if len(line) > 0]

# deprioritize jobs
cmd = 'qalter -p -1023 %s' %(','.join(jobs))
os.system(cmd)

# get interactive node
cmd = 'ish -q interactive -l h_vmem=64g -l h_rt=36000 "qalter -p 0 %s"' %(','.join(jobs))
print cmd

