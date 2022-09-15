__author__ = 'alipirani'
import subprocess
import re
import os
import errno
import glob
from datetime import datetime
from config_settings import ConfigSectionMap
from logging_subprocess import *
from log_modules import keep_logging
from modules.generate_cluster_jobs import *

""" Raw data Pre-processing using Trimmomatic """
def clean_reads(filenames_array, type, out_path, logger, Config, scheduler):
    for file in filenames_array:
        file2 = file.replace('_R1.fastq.gz', '_R2.fastq.gz')
        forward_paired = "/tmp/%s_" % (os.path.basename(file)).replace('_R1.fastq.gz', '') + ConfigSectionMap("Trimmomatic", Config)['f_p']
        reverse_paired = "/tmp/%s_" % (os.path.basename(file)).replace('_R1.fastq.gz', '') + ConfigSectionMap("Trimmomatic", Config)['r_p']
        forward_unpaired = "/tmp/%s_" % (os.path.basename(file)).replace('_R1.fastq.gz', '') + ConfigSectionMap("Trimmomatic", Config)['f_up']
        reverse_unpaired = "/tmp/%s_" % (os.path.basename(file)).replace('_R1.fastq.gz', '') + ConfigSectionMap("Trimmomatic", Config)['r_up']
        adapter_file = ConfigSectionMap("Trimmomatic", Config)['adaptor_filepath']
        clean_filenames = forward_paired + " " + forward_unpaired + " " + reverse_paired + " " + reverse_unpaired
        illumina_string = 'ILLUMINACLIP:' + adapter_file + ConfigSectionMap("Trimmomatic", Config)['colon'] + ConfigSectionMap("Trimmomatic", Config)['seed_mismatches'] + ConfigSectionMap("Trimmomatic", Config)['colon'] + ConfigSectionMap("Trimmomatic", Config)['palindrome_clipthreshold'] + ConfigSectionMap("Trimmomatic", Config)['colon'] + ConfigSectionMap("Trimmomatic", Config)['simple_clipthreshold']
        sliding_string = 'SLIDINGWINDOW:' + ConfigSectionMap("Trimmomatic", Config)['window_size'] + ConfigSectionMap("Trimmomatic", Config)['colon'] + ConfigSectionMap("Trimmomatic", Config)['window_size_quality']
        minlen_string = 'MINLEN:' + ConfigSectionMap("Trimmomatic", Config)['minlength']
        headcrop_string = 'HEADCROP:' + ConfigSectionMap("Trimmomatic", Config)['headcrop_length']


        
        cmdstring = "trimmomatic PE " + file + " " + file2 + " " + clean_filenames + " " + illumina_string + " " + sliding_string + " " + minlen_string + " 2> %s/%s_trim_out.log" % (out_path, os.path.basename(os.path.basename(file)))
        keep_logging(cmdstring, cmdstring, logger, 'info')

        fastqc_forward_cmd = "%s -o %s %s" % (ConfigSectionMap("fastqc", Config)['base_cmd'], out_path, forward_paired)
        fastqc_reverse_cmd = "%s -o %s %s" % (ConfigSectionMap("fastqc", Config)['base_cmd'], out_path, reverse_paired)

        keep_logging(fastqc_forward_cmd, cmdstring, logger, 'info')
        keep_logging(fastqc_reverse_cmd, cmdstring, logger, 'info')

        file_prefix = out_path + "/" + os.path.basename(file.replace('_R1.fastq.gz',''))
        trim_commands = cmdstring + "\n" + fastqc_forward_cmd + "\n" + fastqc_forward_cmd
        job_filename = generate_cluster_jobs(trim_commands, file_prefix, scheduler, Config, logger)
          

        # try:
        #     call(cmdstring, logger)
        #     call(fastqc_forward_cmd, logger)
        #     call(fastqc_reverse_cmd, logger)
        #     #print ""
        # except sp.CalledProcessError:
        #     keep_logging('Error in Trimmomatic Pre-processing step. Exiting.', 'Error in Trimmomatic Pre-processing step. Exiting.', logger, 'exception')
        #     sys.exit(1)
    