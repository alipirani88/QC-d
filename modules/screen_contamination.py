__author__ = 'alipirani'
import os
import subprocess
from modules.log_modules import keep_logging
from modules.logging_subprocess import *
from config_settings import ConfigSectionMap
from modules.generate_cluster_jobs import *

def screen_contamination(filenames_array, Config, logger, output_folder, type, samples, fastq_screen_directory, cluster):
    if type == "PE":
        for file in filenames_array:
            file_prefix = fastq_screen_directory + "/" + os.path.basename(file).replace('.gz', '')
            #fastq_screen_forward_cmd = "%s/%s/%s --subset %s --force --outdir %s --aligner %s %s %s" % (ConfigSectionMap("bin_path", Config)['binbase'], ConfigSectionMap("fastq_screen", Config)['fastq_screen_bin'], ConfigSectionMap("fastq_screen", Config)['base_cmd'], ConfigSectionMap("fastq_screen", Config)['subset'], fastq_screen_directory, ConfigSectionMap("fastq_screen", Config)['aligner'], file, file.replace('_R1_', '_R2_'))
            fastq_screen_forward_cmd = "%s/%s/%s --subset %s --force --outdir %s --aligner %s %s" % (ConfigSectionMap("bin_path", Config)['binbase'], ConfigSectionMap("fastq_screen", Config)['fastq_screen_bin'], ConfigSectionMap("fastq_screen", Config)['base_cmd'], ConfigSectionMap("fastq_screen", Config)['subset'], fastq_screen_directory, ConfigSectionMap("fastq_screen", Config)['aligner'], file)
            keep_logging(fastq_screen_forward_cmd, fastq_screen_forward_cmd, logger, 'debug')
            if cluster == "cluster":
                generate_cluster_jobs(fastq_screen_forward_cmd, file_prefix, Config, logger)
            elif cluster == "local":
                call(fastq_screen_forward_cmd, logger)
    elif type == "SE":
        for file in filenames_array:
            file_prefix = fastq_screen_directory + "/" + os.path.basename(file).replace('.gz', '')
            fastq_screen_forward_cmd = "%s/%s/%s --subset %s --force --outdir %s --aligner %s %s" % (ConfigSectionMap("bin_path", Config)['binbase'], ConfigSectionMap("fastq_screen", Config)['fastq_screen_bin'], ConfigSectionMap("fastq_screen", Config)['base_cmd'], ConfigSectionMap("fastq_screen", Config)['subset'], fastq_screen_directory, ConfigSectionMap("fastq_screen", Config)['aligner'], file)
            keep_logging(fastq_screen_forward_cmd, fastq_screen_forward_cmd, logger, 'debug')
            if cluster == "cluster":
                generate_cluster_jobs(fastq_screen_forward_cmd, file_prefix, Config, logger)
            elif cluster == "local":
                call(fastq_screen_forward_cmd, logger)