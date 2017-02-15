__author__ = 'alipirani'

import os
import readline
import argparse
from joblib import Parallel, delayed
import multiprocessing

def run_parallel(parallel_local_cmds):
    complete = 0
    num_cores = multiprocessing.cpu_count()
    print "The option cluster was set to parallel-local.\n The pipeline will use all the available cores to run the jobs\n"
    print "No of cores available: %s\n" % num_cores
    results = Parallel(n_jobs=num_cores)(delayed(run_command)(i) for i in parallel_local_cmds)
    complete = 1
    return complete

def run_command(i):
    os.system(i)
    done = "Running: %s\n" % i
    print done
    return done