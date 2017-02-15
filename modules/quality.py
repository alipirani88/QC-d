__author__ = 'alipirani'

import os
import subprocess
from modules.log_modules import keep_logging
from modules.logging_subprocess import *
from config_settings import ConfigSectionMap
from itertools import izip

def quality(filenames_array, Config, logger, output_folder, type, samples, fastqc_forward_directory, fastqc_reverse_directory):
    if type == "PE":
        for file in filenames_array:
            fastqc_forward_cmd = "%s/%s/%s -o %s %s" % (ConfigSectionMap("bin_path", Config)['binbase'], ConfigSectionMap("fastqc", Config)['fastqc_bin'], ConfigSectionMap("fastqc", Config)['base_cmd'], fastqc_forward_directory, file)
            fastqc_reverse_cmd = "%s/%s/%s -o %s %s" % (ConfigSectionMap("bin_path", Config)['binbase'], ConfigSectionMap("fastqc", Config)['fastqc_bin'], ConfigSectionMap("fastqc", Config)['base_cmd'], fastqc_reverse_directory, file.replace('_R1_', '_R2_'))
            keep_logging(fastqc_forward_cmd, fastqc_forward_cmd, logger, 'debug')
            keep_logging(fastqc_reverse_cmd, fastqc_reverse_cmd, logger, 'debug')
            call(fastqc_forward_cmd, logger)
            call(fastqc_reverse_cmd, logger)
    elif type == "SE":
        for file in filenames_array:
            fastqc_forward_cmd = "%s/%s/%s -o %s %s" % (ConfigSectionMap("bin_path", Config)['binbase'], ConfigSectionMap("fastqc", Config)['fastqc_bin'], ConfigSectionMap("fastqc", Config)['base_cmd'], fastqc_forward_directory, file)
            keep_logging(fastqc_forward_cmd, fastqc_forward_cmd, logger, 'debug')
            call(fastqc_forward_cmd, logger)

