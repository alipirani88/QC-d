__author__ = 'alipirani'
import os
import subprocess
from modules.log_modules import keep_logging
from modules.logging_subprocess import *
from config_settings import ConfigSectionMap

def generate_cluster_jobs(cmd, jobname, Config, logger):
    job_filename = jobname + ".pbs"
    job = os.path.basename(jobname)
    #cluster_parameters = ConfigSectionMap("cluster", Config)['cluster_parameters']
    #cluster_email = ConfigSectionMap("cluster", Config)['cluster_email']
    cluster_resources = ConfigSectionMap("cluster", Config)['cluster_resources']
    #cluster_account = ConfigSectionMap("cluster", Config)['cluster_account']
    cluster_account = "#PBS -q fluxod\n#PBS -A esnitkin_fluxod\n#PBS -l qos=flux"
    cluster_email = "#PBS -M apirani@med.umich.edu"
    cluster_parameters = "#PBS -m abe\n#PBS -V"
    job_message = "Generating PBS Scripts for Job: %s" % job
    keep_logging(job_message, job_message, logger, 'info')
    with open(job_filename, 'w') as out:
        job_title = "#PBS -N %s" % job
        out.write(job_title+'\n')
        out.write(cluster_parameters+'\n')
        out.write(cluster_email+'\n')
        out.write(cluster_resources+'\n')
        out.write(cluster_account+'\n')
        out.write('\n'+cmd+'\n')