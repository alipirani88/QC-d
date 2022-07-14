__author__ = 'alipirani'
import os
import subprocess
#import statistics
from modules.log_modules import keep_logging
from logging_subprocess import *
from config_settings import ConfigSectionMap
from itertools import izip
from modules.generate_cluster_jobs import *

def multiqc(analysis_folder, filename, Config, logger, Multiqc_reports_directory, cluster, scheduler):
    message = "Running MultiQC on %s" % analysis_folder
    keep_logging('', message, logger, 'debug')
    run_multiqc_cmd = "%s %s --force --filename %s --outdir %s" % (ConfigSectionMap("multiqc", Config)['base_cmd'], analysis_folder, filename, Multiqc_reports_directory)
    keep_logging(run_multiqc_cmd, run_multiqc_cmd, logger, 'debug')
    
    if cluster == "cluster":
                jobpath = "%s/multiqc.sbat" % analysis_folder
                job_filename = generate_cluster_jobs(run_multiqc_cmd, jobpath, scheduler, Config, logger)
                keep_logging('Generating Job - %s' % job_filename, 'Generating Job - %s' % job_filename, logger, 'debug')
                #call("sbatch %s" % job_filename, logger)
                
    elif cluster == "local":
        call(run_multiqc_cmd, logger)
                