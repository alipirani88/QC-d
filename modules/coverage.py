__author__ = 'alipirani'
import os
import subprocess
import statistics
from modules.log_modules import keep_logging
from logging_subprocess import *
from config_settings import ConfigSectionMap
from itertools import izip

def coverage(filenames_array, Config, logger, output_folder, type, samples, size, prefix):
    temp_forward_coverage = output_folder + "/%s_temp_forward.txt" % prefix
    temp_reverse_coverage = output_folder + "/%s_temp_reverse.txt" % prefix
    f1=open(temp_forward_coverage, 'w+')
    f2=open(temp_reverse_coverage, 'w+')
    if type == "PE":
        for file in filenames_array:
            if file.endswith('.gz'):
                coverage_msg_forward = "Calculating coverage for file: %s\n" % file
                coverage_msg_reverse = "Calculating coverage for file: %s\n" % file.replace('_R1_', '_R2_')
                keep_logging('', coverage_msg_forward, logger, 'debug')
                keep_logging('', coverage_msg_reverse, logger, 'debug')
                coverage_cmd_forward = "zcat %s | %s" % (file, ConfigSectionMap("coverage", Config)['awk_cmd'])
                coverage_cmd_reverse = "zcat %s | %s" % (file.replace('_R1_', '_R2_'), ConfigSectionMap("coverage", Config)['awk_cmd'])
            else:
                coverage_msg_forward = "Calculating coverage for file: %s\n" % file
                coverage_msg_reverse = "Calculating coverage for file: %s\n" % file.replace('_R1_', '_R2_')
                keep_logging('', coverage_msg_forward, logger, 'debug')
                keep_logging('', coverage_msg_reverse, logger, 'debug')
                coverage_cmd_forward = "cat %s | %s" % (file, ConfigSectionMap("coverage", Config)['awk_cmd'])
                coverage_cmd_reverse = "cat %s | %s" % (file.replace('_R1_', '_R2_'), ConfigSectionMap("coverage", Config)['awk_cmd'])
            keep_logging('', coverage_cmd_forward, logger, 'debug')
            keep_logging('', coverage_cmd_reverse, logger, 'debug')
            proc = subprocess.Popen([coverage_cmd_forward], stdout=subprocess.PIPE, shell=True)
            (out2, err2) = proc.communicate()
            f1.write(out2)
            proc2 = subprocess.Popen([coverage_cmd_reverse], stdout=subprocess.PIPE, shell=True)
            (out, err) = proc2.communicate()
            f2.write(out)
        temp_final_file = "%s/%s_temp_final.txt" % (output_folder, prefix)
        sample_file = "%s/%s" % (output_folder, os.path.basename(samples))
        with open(temp_final_file, 'w') as res, open(sample_file) as f1, open(temp_forward_coverage) as f2, open(temp_reverse_coverage) as f3:
            for line1, line2, line3 in zip(f1, f2, f3):
                res.write("{}\t{}\t{}\n".format(line1.rstrip(), line2.rstrip(), line3.rstrip()))
        final_coverage_file = "%s/%s_Final_Coverage.txt" % (output_folder, prefix)
        f3=open(final_coverage_file, 'w+')
        header = "Sample_name\tForward_reads\tForward_unique_reads\tPerc_Forward_unique_reads\tMost_abundant_Forward_read(adapter)\tMost_abundant_Forward_read_count\tPerc_Most_abundant_Forward_read\tAvg_Forward_readLength\tReverse_reads\tReverse_unique_reads\tPerc_Reverse_unique_reads\tMost_abundant_Reverse_read(adapter)\tMost_abundant_Reverse_read_count\tPerc_Most_abundant_Reverse_read\tAvg_Reverse_readLength\tCoverage\n"
        f3.write(header)
        with open(temp_final_file, 'a+') as fp:
            for line in fp:
                line = line.strip()
                line_split = line.split('\t')
                avg_read_length = statistics.mean([float(line_split[7]), float(line_split[14])])
                final_coverage = (int(line_split[1]) * 2 * avg_read_length) / int(size)
                #print final_coverage
                print_string = line + "\t" + str(final_coverage) + "\n"
                f3.write(print_string)
        f3.close()

    elif type == "SE":
        for file in filenames_array:
            coverage_msg_forward = "Calculating coverage for file: %s\n" % file
            if file.endswith('.gz'):
                coverage_msg_forward = "Calculating coverage for file: %s\n" % file
                keep_logging('', coverage_msg_forward, logger, 'debug')
                coverage_cmd_forward = "zcat %s | %s" % (file, ConfigSectionMap("coverage", Config)['awk_cmd'])
            else:
                coverage_msg_forward = "Calculating coverage for file: %s\n" % file
                keep_logging('', coverage_msg_forward, logger, 'debug')
                coverage_cmd_forward = "cat %s | %s" % (file, ConfigSectionMap("coverage", Config)['awk_cmd'])
            keep_logging('', coverage_cmd_forward, logger, 'debug')
            proc = subprocess.Popen([coverage_cmd_forward], stdout=subprocess.PIPE, shell=True)
            (out2, err2) = proc.communicate()
            f1.write(out2)
        temp_final_file = "%s/%s_temp_final.txt" % (output_folder, prefix)
        sample_file = "%s/%s" % (output_folder, os.path.basename(samples))
        with open(temp_final_file, 'w') as res, open(sample_file) as f1, open(temp_forward_coverage) as f2:
            for line1, line2 in zip(f1, f2):
                res.write("{}\t{}\n".format(line1.rstrip(), line2.rstrip()))
        final_coverage_file = "%s/%s_Final_Coverage.txt" % (output_folder, prefix)
        f3=open(final_coverage_file, 'w+')
        header = "Sample_name\tForward_reads\tForward_unique_reads\tPerc_Forward_unique_reads\tMost_abundant_Forward_read(adapter)\tMost_abundant_Forward_read_count\tPerc_Most_abundant_Forward_read\tAvg_Forward_readLength\tCoverage\n"
        f3.write(header)
        with open(temp_final_file, 'a+') as fp:
            for line in fp:
                line = line.strip()
                line_split = line.split('\t')
                avg_read_length = float(line_split[7])
                final_coverage = (int(line_split[1]) * avg_read_length) / int(size)
                #print final_coverage
                print_string = line + "\t" + str(final_coverage) + "\n"
                f3.write(print_string)
        f3.close()

    os.system("rm %s %s %s" % (temp_forward_coverage, temp_reverse_coverage, temp_final_file))
    keep_logging('Coverage Report - %s\n' % final_coverage_file, 'Coverage Report - %s\n' % final_coverage_file, logger, 'info')