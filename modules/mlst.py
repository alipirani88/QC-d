__author__ = 'alipirani'

import os
import subprocess
import re
from modules.log_modules import keep_logging
from logging_subprocess import *
from config_settings import ConfigSectionMap
from itertools import izip
from modules.generate_cluster_jobs import *


def mlst(filenames_array, Config, logger, output_folder, type, samples, mlst_directory, cluster, scheduler, mlstdb):
    if type == "PE":
        for file in filenames_array:
            if re.search('R1_001_final.fastq.gz', file):
                second_part = file.replace("R1_001_final.fastq.gz", "R2_001_final.fastq.gz")
                first_part_split = file.split('R1_001_final.fastq.gz')
                first_part = first_part_split[0].replace('_L001', '')
                first_part = re.sub("_S[0-9].*_", "", first_part)
            elif re.search('_R1.fastq.gz', file):
                second_part = file.replace("_R1.fastq.gz", "_R2.fastq.gz")
                first_part_split = file.split('_R1.fastq.gz')
                first_part = first_part_split[0].replace('_L001', '')
                first_part = re.sub("_S[0-9].*_", "", first_part)
                # Changed on 03/15/2019
            elif re.search('R1.fastq.gz', file):
                second_part = file.replace("R1.fastq.gz", "R2.fastq.gz")
                first_part_split = file.split('R1.fastq.gz')
                first_part = first_part_split[0].replace('_L001', '')
                first_part = re.sub("_S.*_", "", first_part)
                # Changed on 03/15/2019
                first_part = re.sub("_S[0-9].*_", "", first_part)
            elif re.search('1_combine.fastq.gz', file):
                second_part = file.replace("1_combine.fastq.gz", "2_combine.fastq.gz")
                first_part_split = file.split('1_combine.fastq.gz')
                first_part = first_part_split[0].replace('_L001', '')
                first_part = re.sub("_S[0-9].*_", "", first_part)
            elif re.search('1_sequence.fastq.gz', file):
                second_part = file.replace("1_sequence.fastq.gz", "2_sequence.fastq.gz")
                first_part_split = file.split('1_sequence.fastq.gz')
                first_part = first_part_split[0].replace('_L001', '')
                first_part = re.sub("_S[0-9].*_", "", first_part)
            elif re.search('_forward.fastq.gz', file):
                second_part = file.replace("_forward.fastq.gz", "_reverse.fastq.gz")
                first_part_split = file.split('_forward.fastq.gz')
                first_part = first_part_split[0].replace('_L001', '')
                first_part = re.sub("_S[0-9].*_", "", first_part)
            elif re.search('R1_001.fastq.gz', file):
                second_part = file.replace("R1_001.fastq.gz", "R2_001.fastq.gz")
                first_part_split = file.split('R1_001.fastq.gz')
                first_part = first_part_split[0].replace('_L001', '')
                first_part = re.sub("_S[0-9].*_", "", first_part)
            elif re.search('_1.fastq.gz', file):
                second_part = file.replace("_1.fastq.gz", "_2.fastq.gz")
                first_part_split = file.split('_1.fastq.gz')
                first_part = first_part_split[0].replace('_L001', '')
                first_part = re.sub("_S[0-9].*_", "", first_part)
            elif re.search('.1.fastq.gz', file):
                second_part = file.replace(".1.fastq.gz", ".2.fastq.gz")
                first_part_split = file.split('.1.fastq.gz')
                first_part = re.sub("_S[0-9].*_", "", first_part)
                #first_part = re.sub("_S.*_", "", first_part)
            # mlst_cmd = "ariba run --verbose --force %s %s %s %s/%s --tmp_dir /tmp/" % (ConfigSectionMap("ariba", Config)['mlst_db_path'], file, file.replace('_R1_', '_R2_'), output_folder, os.path.basename(file)[0:20])
            mlst_cmd = "ariba run --verbose --force %s %s %s %s/%s --tmp_dir /tmp/" % (mlstdb, file, second_part, output_folder, os.path.basename(first_part))
            #keep_logging(mlst_cmd, mlst_cmd, logger, 'debug')
            if cluster == "cluster":
                job_prefix = "%s/%s" % (output_folder, os.path.basename(first_part))
                generate_cluster_jobs(mlst_cmd, job_prefix, scheduler, Config, logger)

    elif type == "SE":
        keep_logging('Ariba requires PE files', 'Ariba requires PE files', logger, 'debug')
