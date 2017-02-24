__author__ = 'alipirani'
import os
import subprocess
import statistics
from modules.log_modules import keep_logging
from logging_subprocess import *
from config_settings import ConfigSectionMap
from itertools import izip

def multiqc(analysis_folder, filename, Config, logger, Multiqc_reports_directory):
    message = "Running MultiQC on %s" % analysis_folder
    keep_logging(message, message, logger, 'info')
    run_multiqc_cmd = "%s/%s %s --force --filename %s --outdir %s" % (ConfigSectionMap("multiqc", Config)['multiqc_bin'], ConfigSectionMap("multiqc", Config)['base_cmd'], analysis_folder, filename, Multiqc_reports_directory)
    keep_logging('', run_multiqc_cmd, logger, 'debug')
    call(run_multiqc_cmd, logger)
