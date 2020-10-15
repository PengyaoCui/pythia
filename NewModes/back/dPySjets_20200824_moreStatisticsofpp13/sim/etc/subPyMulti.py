#!/usr/bin/env python

#SBATCH -J PyMulti
#SBATCH -o log/%x_%A_%a.out
#SBATCH -e log/%x_%A_%a.err
#SBATCH --array=1-50
#SBATCH --ntasks=4
###############################################################################

import sys, os
import subprocess
import multiprocessing as mp
###############################################################################

sys.path.append(os.getcwd())
###############################################################################

mod_inst = 'local/opt/modules/master'
mod_path = os.environ['HOME'] + '/' + mod_inst
mod_init = '{}/init/python.py'.format(mod_path)
exec(open(mod_init).read())
###############################################################################

exef = 'PyMulti'
###############################################################################

def load_env_mod() :
    if not module('is-loaded', 'batch') :
        module('use', 'etc/modulefiles')
        module('load', 'batch')
###############################################################################

def run_exec(tune) :
    run = subprocess.run([exef,
        '{}'.format(tune),
        '{}'.format(os.environ['SLURM_ARRAY_JOB_ID']),
        '{}'.format(os.environ['SLURM_ARRAY_TASK_ID']),
         os.environ['SLURM_JOB_NAME'] ])
###############################################################################

if __name__ == '__main__' :
    load_env_mod()

    tunes = [ i for i in range(4) ]
    with mp.Pool(len(tunes)) as p :
        p.map(run_exec, tunes)
###############################################################################
