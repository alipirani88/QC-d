__author__ = 'alipirani'
import os
import subprocess
from modules.log_modules import keep_logging
from logging_subprocess import *
from config_settings import ConfigSectionMap
from modules.generate_cluster_jobs import *
from modules.run_parallel import *


def kraken_report(filenames_array, Config, logger, output_folder, type, samples, kraken_directory, cluster, scheduler):
    kraken_report_array = []
    kraken_report = open("%s/Kraken_report.sh" % (kraken_directory), 'w+')
    echo = "echo \"Sample,Percentage of reads for Species,# of reads for Species, Species\" > %s/Kraken_report_final.csv" % (kraken_directory)
    
    prepare_report = "for i in %s/*_report.txt; do grep -w \'S\' $i | sort -k1n | tail -n1; done > %s/Kraken_report_temp.txt\npaste %s %s/Kraken_report_temp.txt > %s/Kraken_report_combined.txt\n" \
                             "awk -F\'\\t\' \'BEGIN{OFS=\",\"};{print $1, $2, $3, $7}\' %s/Kraken_report_combined.txt >> %s/Kraken_report_final.csv\n" \
                             "sed -i \'s/\\s//g\' %s/Kraken_report_final.csv" % (kraken_directory, kraken_directory, samples, kraken_directory, kraken_directory, kraken_directory, kraken_directory, kraken_directory)
    subprocess.call(echo, shell=True)
    subprocess.call(prepare_report, shell=True)

    kraken_report.write("%s\n%s" % (echo, prepare_report))
    keep_logging(echo, echo, logger, 'debug')
    keep_logging(prepare_report, prepare_report, logger, 'debug')
    
    keep_logging("Writing Kraken Report commands to - %s/Kraken_report.sh" % kraken_directory, "Writing Kraken Report commands to - %s/Kraken_report.sh" % kraken_directory, logger, 'debug')
    
    keep_logging('',
                 "\nKraken Report generated with Kraken_report.sh - %s/Kraken_report_final.csv" % kraken_directory, logger, 'debug')
