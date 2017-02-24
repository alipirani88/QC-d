__author__ = 'alipirani'
import os
import subprocess
from modules.log_modules import keep_logging
from logging_subprocess import *
from config_settings import ConfigSectionMap
from modules.generate_cluster_jobs import *
from modules.run_parallel import *


def kraken_report(filenames_array, Config, logger, output_folder, type, samples, kraken_directory, cluster):
    kraken_report_array = []
    echo = "echo \"Sample,Percentage of reads for Species,# of reads for Species, Species\" > %s/Kraken_report_final.csv" % (kraken_directory)
    os.system(echo)
    for file in filenames_array:
        file_prefix = kraken_directory + "/" + os.path.basename(file)[0:20]
        kraken_out = file_prefix + "_kraken"
        report_cmd = "%s/%s/kraken-report --db %s %s > %s_report.txt" % (ConfigSectionMap("bin_path", Config)['binbase'], ConfigSectionMap("kraken", Config)['kraken_bin'], ConfigSectionMap("kraken", Config)['db_path'], kraken_out, kraken_out)
        keep_logging('', report_cmd, logger, 'debug')
        if cluster == "cluster":
            generate_cluster_jobs(report_cmd, file_prefix, Config, logger)
        elif cluster == "parallel-local":
            kraken_report_array.append(report_cmd)
        elif cluster == "local":
            call(report_cmd)
    if cluster == "parallel-local":
        complete = run_parallel(kraken_report_array)
        prepare_report = "for i in %s/*_report.txt; do grep -w \'S\' $i | sort -k1n | tail -n1; done > %s/Kraken_report.txt\npaste %s %s/Kraken_report.txt > %s/Kraken_report_combined.txt\n" \
                             "awk -F\'\\t\' \'BEGIN{OFS=\",\"};{print $1, $2, $3, $7}\' %s/Kraken_report_combined.txt >> %s/Kraken_report_final.csv" % (kraken_directory, kraken_directory, samples, kraken_directory, kraken_directory, kraken_directory, kraken_directory)

        subprocess.call(["for i in %s/*_report.txt; do grep -w 'S' $i | sort -k1n | tail -n1; done > %s/Kraken_report.txt" % (kraken_directory, kraken_directory)], shell=True)
        subprocess.call(["paste %s %s/Kraken_report.txt > %s/Kraken_report_combined.txt" % (samples, kraken_directory, kraken_directory)], shell=True)
        subprocess.call(["awk -F'\t' 'BEGIN{OFS=\",\"};{print $1, $2, $3, $7}' %s/Kraken_report_combined.txt >> %s/Kraken_report_final.csv" % (kraken_directory, kraken_directory)], shell=True)

        #print "Running:\n%s" % prepare_report
        keep_logging('', prepare_report, logger, 'debug')

    prepare_report = "for i in %s/*_report.txt; do grep -w \'S\' $i | sort -k1n | tail -n1; done > %s/Kraken_report.txt\npaste %s %s/Kraken_report.txt > %s/Kraken_report_combined.txt\n" \
                             "awk -F\'\\t\' \'BEGIN{OFS=\",\"};{print $1, $2, $3, $7}\' %s/Kraken_report_combined.txt >> %s/Kraken_report_final.csv\n" \
                             "sed -i \'s/\\s//g\' %s/Kraken_report_final.csv" % (kraken_directory, kraken_directory, samples, kraken_directory, kraken_directory, kraken_directory, kraken_directory, kraken_directory)
    #print "Running:\n%s" % prepare_report
    keep_logging('', prepare_report, logger, 'debug')