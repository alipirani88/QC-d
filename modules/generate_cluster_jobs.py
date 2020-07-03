__author__ = 'alipirani'
import os
import subprocess
from modules.log_modules import keep_logging
from logging_subprocess import *
from config_settings import ConfigSectionMap



def get_scheduler_directive(scheduler, Config):
    """ Generate Cluster Directive lines for a scheduler provided with args.scheduler"""
    # Scheduler Changes here; current changes
    if scheduler and scheduler == "SLURM":
        script_Directive = "#SBATCH"
        job_name_flag = "--job-name="
        scheduler_directives = "#SBATCH --mail-user=%s\n#SBATCH --mail-type=%s\n#SBATCH --export=ALL\n#SBATCH --partition=%s\n#SBATCH --account=%s\n#SBATCH %s\n" \
                          % (ConfigSectionMap("slurm", Config)['email'],
                             ConfigSectionMap("slurm", Config)['notification'],
                             ConfigSectionMap("slurm", Config)['partition'],
                             ConfigSectionMap("slurm", Config)['flux_account'],
                             ConfigSectionMap("slurm", Config)['resources'])
    elif scheduler and scheduler == "PBS":
        script_Directive = "#PBS"
        job_name_flag = "-N"
        scheduler_directives = "#PBS -M %s\n#PBS -m %s\n#PBS -V\n#PBS -l %s\n#PBS -q %s\n#PBS -A %s\n#PBS -l qos=flux\n" \
                          % (ConfigSectionMap("scheduler", Config)['email'],
                             ConfigSectionMap("scheduler", Config)['notification'],
                             ConfigSectionMap("scheduler", Config)['resources'],
                             ConfigSectionMap("scheduler", Config)['queue'],
                             ConfigSectionMap("scheduler", Config)['flux_account'])
    else:
        script_Directive = "#SBATCH"
        job_name_flag = "--job-name="
        scheduler_directives = "#SBATCH --mail-user=%s\n#SBATCH --mail-type=%s\n#SBATCH --export=ALL\n#SBATCH --partition=%s\n#SBATCH --account=%s\n#SBATCH %s\n" \
                               % (ConfigSectionMap("slurm", Config)['email'],
                                  ConfigSectionMap("slurm", Config)['notification'],
                                  ConfigSectionMap("slurm", Config)['partition'],
                                  ConfigSectionMap("slurm", Config)['flux_account'],
                                  ConfigSectionMap("slurm", Config)['resources'])
    return scheduler_directives, script_Directive, job_name_flag

def generate_cluster_jobs(cmd, jobname, scheduler, Config, logger):
    print "Generating Cluster jobs for %s\n" % jobname
    scheduler_directives, script_Directive, job_name_flag = get_scheduler_directive(scheduler, Config)
    if scheduler == "SLURM":
        job_filename = jobname + ".sbat"
    elif scheduler == "PBS":
        job_filename = jobname + ".pbs"
    else:
        print "Provide Scheduler to generate cluster jobs"
        exit()
    with open(job_filename, 'w') as out:
        job_title = "%s %s%s" % (script_Directive, job_name_flag, jobname)
        out.write("#!/bin/sh" + '\n')
        out.write(script_Directive + " " + job_name_flag + os.path.basename(job_title) + '\n')
        out.write(scheduler_directives + '\n')
        out.write(cmd + '\n')

    return job_filename
