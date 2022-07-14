__author__ = 'alipirani'

import os
import subprocess
from modules.log_modules import keep_logging
from logging_subprocess import *
from config_settings import ConfigSectionMap
from itertools import izip
from modules.generate_cluster_jobs import *

def quality(filenames_array, Config, logger, output_folder, type, samples, fastqc_main_directory, fastqc_forward_directory, fastqc_reverse_directory, cluster, scheduler):
    if type == "PE":
        fastqc_fwd = open("%s/fastqc.sh" % (fastqc_main_directory), 'w+')
        for file in filenames_array:
            fastqc_msg_forward = "Running FastQC on Forward-end file: %s\n" % file
            fastqc_msg_reverse = "Running FastQC on Reverse-end file: %s\n" % file.replace('_R1_', '_R2_')
            keep_logging('', fastqc_msg_forward, logger, 'debug')
            keep_logging('', fastqc_msg_reverse, logger, 'debug')
            
            fastqc_forward_cmd = "%s -o %s %s" % (ConfigSectionMap("fastqc", Config)['base_cmd'], fastqc_forward_directory, file)
            fastqc_reverse_cmd = "%s -o %s %s" % (ConfigSectionMap("fastqc", Config)['base_cmd'], fastqc_reverse_directory, file.replace('_R1', '_R2'))
            keep_logging('', fastqc_forward_cmd, logger, 'debug')
            keep_logging('', fastqc_reverse_cmd, logger, 'debug')
            fastqc_fwd.write(fastqc_forward_cmd + ' && ' + fastqc_reverse_cmd + '\n')
            
            if cluster == "cluster":
                jobpath = "%s/%s" % (fastqc_main_directory, os.path.basename(file.replace('_R1_001.fastq.gz', '')))
                job_filename = generate_cluster_jobs(fastqc_forward_cmd + ' && ' + fastqc_reverse_cmd, jobpath, scheduler, Config, logger)
                keep_logging('Generating Job - %s' % job_filename, 'Generating Job - %s' % job_filename, logger, 'debug')
                #call("sbatch %s" % job_filename, logger)
                
            elif cluster == "local":
                call(fastqc_forward_cmd, logger)
                call(fastqc_reverse_cmd, logger)
        
    elif type == "SE":
        for file in filenames_array:
            fastqc_msg_forward = "Running FastQC on Forward-end file: %s\n" % file
            keep_logging('', fastqc_msg_forward, logger, 'debug')
            fastqc_forward_cmd = "%s -o %s %s" % (ConfigSectionMap("fastqc", Config)['base_cmd'], fastqc_forward_directory, file)
            keep_logging('', fastqc_forward_cmd, logger, 'debug')
            call(fastqc_forward_cmd, logger)
