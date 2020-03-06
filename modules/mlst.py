__author__ = 'alipirani'

import os
import subprocess
from modules.log_modules import keep_logging
from logging_subprocess import *
from config_settings import ConfigSectionMap
from itertools import izip
from modules.generate_cluster_jobs import *

def mlst(filenames_array, Config, logger, output_folder, type, samples, mlst_directory, cluster, scheduler):
    if type == "PE":
        for file in filenames_array:
            mlst_cmd = "ariba run --verbose --force %s %s %s %s/%s --tmp_dir /tmp/" % (ConfigSectionMap("ariba", Config)['mlst_db_path'], file, file.replace('_R1_', '_R2_'), output_folder, os.path.basename(file)[0:20])
            #keep_logging(mlst_cmd, mlst_cmd, logger, 'debug')
            if cluster == "cluster":
                job_prefix = "%s/%s" % (output_folder, os.path.basename(file)[0:20])
                generate_cluster_jobs(mlst_cmd, job_prefix, scheduler, Config, logger)

    elif type == "SE":
        keep_logging('Ariba requires PE files', 'Ariba requires PE files', logger, 'debug')
