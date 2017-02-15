__author__ = 'alipirani'
import os
import subprocess
from modules.log_modules import keep_logging
from modules.logging_subprocess import *
from config_settings import ConfigSectionMap
from modules.generate_cluster_jobs import *
from modules.bwa import align_bwa
import gzip
import re
from modules.samtools import *
from modules.gatk import *


def coverage_depth_analysis(filenames_array, Config, logger, output_folder, type, samples, coverage_depth_directory, cluster, reference):
    files_to_delete = []
    command_list = []
    if type == "PE":
        for file in filenames_array:
            filename_base = os.path.basename(file)
            if "R1_001_final.fastq.gz" in filename_base:
                reverse_file = file.replace("R1_001_final.fastq.gz", "R2_001_final.fastq.gz")
                first_part_split = filename_base.split('R1_001_final.fastq.gz')
                first_part = first_part_split[0].replace('_L001', '')
            elif "R1.fastq.gz" in filename_base:
                reverse_file = file.replace("R1.fastq.gz", "R2.fastq.gz")
                first_part_split = filename_base.split('R1.fastq.gz')
                first_part = first_part_split[0].replace('_L001', '')
            elif "1_combine.fastq.gz" in filename_base:
                reverse_file = file.replace("1_combine.fastq.gz", "2_combine.fastq.gz")
                first_part_split = filename_base.split('1_combine.fastq.gz')
                first_part = first_part_split[0].replace('_L001', '')
            elif "1_sequence.fastq.gz" in filename_base:
                reverse_file = file.replace("1_sequence.fastq.gz", "2_sequence.fastq.gz")
                first_part_split = filename_base.split('1_sequence.fastq.gz')
                first_part = first_part_split[0].replace('_L001', '')
            elif "_forward.fastq.gz" in filename_base:
                reverse_file = file.replace("_forward.fastq.gz", "_reverse.fastq.gz")
                first_part_split = filename_base.split('_forward.fastq.gz')
                first_part = first_part_split[0].replace('_L001', '')
            elif "R1_001.fastq.gz" in filename_base:
                reverse_file = file.replace("R1_001.fastq.gz", "R2_001.fastq.gz")
                first_part_split = filename_base.split('R1_001.fastq.gz')
                first_part = first_part_split[0].replace('_L001', '')
            elif "_1.fastq.gz" in filename_base:
                reverse_file = file.replace("_1.fastq.gz", "_2.fastq.gz")
                first_part_split = filename_base.split('_1.fastq.gz')
                first_part = first_part_split[0].replace('_L001', '')
            else:
                print "Using Standard second file naming convention"
                reverse_file = file.replace("_R1_", "_R2_")
                first_part_split = filename_base.split('_R1.fastq.gz')
                first_part = first_part_split[0].replace('_L001', '')
            file_prefix = coverage_depth_directory + "/" + first_part
            analysis = first_part
            if file.endswith('.gz'):
                keep_logging("Generating command list to create cluster jobs", "Generating command list to create cluster jobs", logger, 'info')
                split_field = prepare_readgroup(file, logger)
                command_list, files_to_delete = align_bwa(ConfigSectionMap("bin_path", Config)['binbase'] + ConfigSectionMap("bwa", Config)['bwa_bin'] + ConfigSectionMap("bwa", Config)['base_cmd'],file, reverse_file, coverage_depth_directory, reference, split_field, first_part, files_to_delete, logger, Config, type, command_list)
                out_sam = files_to_delete[0]
                command_list, files_to_delete = samtobam(out_sam, coverage_depth_directory, analysis, files_to_delete, logger, Config, command_list)
                out_bam = files_to_delete[1]
                command_list, files_to_delete = sort_bam(out_bam, coverage_depth_directory, analysis, logger, Config, command_list, files_to_delete)
                out_sort_bam = files_to_delete[2]
                command_list = index_bam(out_sort_bam, coverage_depth_directory, logger, Config, command_list, files_to_delete)
                command_list, gatk_depth_of_coverage_file = gatk_DepthOfCoverage(out_sort_bam, coverage_depth_directory, analysis, reference, logger, Config, command_list)
                command_list = flagstat(out_sort_bam, output_folder, analysis, logger, Config, command_list)
                coverage_depth_cmd = ""

                for i in command_list:
                    coverage_depth_cmd = coverage_depth_cmd + i + "\n"
                keep_logging("The coverage Depth commands for file %s are:\n", "The coverage Depth commands for file %s are:\n", logger, 'info')
                keep_logging(coverage_depth_cmd, coverage_depth_cmd, logger, 'debug')


                if cluster == "cluster":
                    generate_cluster_jobs(coverage_depth_cmd, file_prefix, Config, logger)
                else:
                    f3=open(file_prefix + '_commands.sh', 'w+')
                    f3.write(coverage_depth_cmd)

            else:
                keep_logging("Generating command list to create cluster jobs", "Generating command list to create cluster jobs", logger, 'info')
                split_field = prepare_readgroup(file, logger)
                command_list, files_to_delete = align_bwa(ConfigSectionMap("bin_path", Config)['binbase'] + ConfigSectionMap("bwa", Config)['bwa_bin'],file, reverse_file, coverage_depth_directory, reference, split_field, first_part, files_to_delete, logger, Config, type, command_list)
                out_sam = files_to_delete[0]
                command_list, files_to_delete = samtobam(out_sam, coverage_depth_directory, analysis, files_to_delete, logger, Config, command_list, files_to_delete)
                out_bam = files_to_delete[1]
                command_list, files_to_delete = sort_bam(out_bam, coverage_depth_directory, analysis, logger, Config, command_list, files_to_delete)
                out_sort_bam = files_to_delete[2]
                command_list = index_bam(out_sort_bam, coverage_depth_directory, logger, Config, command_list, files_to_delete)
                command_list, gatk_depth_of_coverage_file = gatk_DepthOfCoverage(out_sorted_bam, coverage_depth_directory, analysis, reference, logger, Config, command_list)
                command_list = flagstat(out_sort_bam, coverage_depth_directory, analysis, logger, Config, command_list)
                coverage_depth_cmd = ""
                for i in command_list:
                    coverage_depth_cmd = coverage_depth_cmd + i + "\n"
                keep_logging("The coverage Depth commands for file %s are:\n", "The coverage Depth commands for file %s are:\n", logger, 'info')
                keep_logging(coverage_depth_cmd, coverage_depth_cmd, logger, 'debug')


                if cluster == "cluster":
                    generate_cluster_jobs(coverage_depth_cmd, file_prefix, Config, logger)
                else:
                    f3=open(file_prefix + '_commands.sh', 'w+')
                    f3.write(coverage_depth_cmd)
    elif type == "SE":
        ###Pending Changes
        for file in filenames_array:
            filename_base = os.path.basename(file)
            if "R1_001_final.fastq.gz" in filename_base:
                reverse_file = file.replace("R1_001_final.fastq.gz", "R2_001_final.fastq.gz")
                first_part_split = filename_base.split('R1_001_final.fastq.gz')
                first_part = first_part_split[0].replace('_L001', '')
            elif "R1.fastq.gz" in filename_base:
                reverse_file = file.replace("R1.fastq.gz", "R2.fastq.gz")
                first_part_split = filename_base.split('R1.fastq.gz')
                first_part = first_part_split[0].replace('_L001', '')
            elif "1_combine.fastq.gz" in filename_base:
                reverse_file = file.replace("1_combine.fastq.gz", "2_combine.fastq.gz")
                first_part_split = filename_base.split('1_combine.fastq.gz')
                first_part = first_part_split[0].replace('_L001', '')
            elif "1_sequence.fastq.gz" in filename_base:
                reverse_file = file.replace("1_sequence.fastq.gz", "2_sequence.fastq.gz")
                first_part_split = filename_base.split('1_sequence.fastq.gz')
                first_part = first_part_split[0].replace('_L001', '')
            elif "_forward.fastq.gz" in filename_base:
                reverse_file = file.replace("_forward.fastq.gz", "_reverse.fastq.gz")
                first_part_split = filename_base.split('_forward.fastq.gz')
                first_part = first_part_split[0].replace('_L001', '')
            elif "R1_001.fastq.gz" in filename_base:
                reverse_file = file.replace("R1_001.fastq.gz", "R2_001.fastq.gz")
                first_part_split = filename_base.split('R1_001.fastq.gz')
                first_part = first_part_split[0].replace('_L001', '')
            elif "_1.fastq.gz" in filename_base:
                reverse_file = file.replace("_1.fastq.gz", "_2.fastq.gz")
                first_part_split = filename_base.split('_1.fastq.gz')
                first_part = first_part_split[0].replace('_L001', '')
            else:
                print "Using Standard second file naming convention"
                reverse_file = file.replace("_R1_", "_R2_")
                first_part_split = filename_base.split('_R1.fastq.gz')
                first_part = first_part_split[0].replace('_L001', '')
            file_prefix = coverage_depth_directory + "/" + first_part

            if file.endswith('.gz'):
                keep_logging("Generating command list to create cluster jobs", "Generating command list to create cluster jobs", logger, 'info')
                split_field = prepare_readgroup(file, logger)
                command_list, files_to_delete = align_bwa(ConfigSectionMap("bin_path", Config)['binbase'] + ConfigSectionMap("bwa", Config)['bwa_bin'],file, reverse_file, coverage_depth_directory, reference, split_field, first_part, files_to_delete, logger, Config, type, command_list)
                out_sam = files_to_delete[0]
                command_list, files_to_delete = samtobam(out_sam, coverage_depth_directory, analysis, files_to_delete, logger, Config, command_list, files_to_delete)
                out_bam = files_to_delete[1]
                command_list, files_to_delete = sort_bam(out_bam, coverage_depth_directory, analysis, logger, Config, command_list, files_to_delete)
                out_sort_bam = files_to_delete[2]
                command_list = index_bam(out_sort_bam, coverage_depth_directory, logger, Config, command_list, files_to_delete)
                command_list, gatk_depth_of_coverage_file = gatk_DepthOfCoverage(out_sorted_bam, coverage_depth_directory, analysis, reference, logger, Config, command_list)
                command_list = flagstat(out_sorted_bam, coverage_depth_directory, analysis, logger, Config, command_list)
                coverage_depth_cmd = ""
                for i in command_list:
                    coverage_depth_cmd = coverage_depth_cmd + i + "\n"
                keep_logging("The coverage Depth commands for file %s are:\n", "The coverage Depth commands for file %s are:\n", logger, 'info')
                keep_logging(coverage_depth_cmd, coverage_depth_cmd, logger, 'debug')


                if cluster == "cluster":
                    generate_cluster_jobs(coverage_depth_cmd, file_prefix, Config, logger)
                else:
                    f3=open(file_prefix + '_commands.sh', 'w+')
                    f3.write(coverage_depth_cmd)

            else:
                keep_logging("Generating command list to create cluster jobs", "Generating command list to create cluster jobs", logger, 'info')
                split_field = prepare_readgroup(file, logger)
                command_list, files_to_delete = align_bwa(ConfigSectionMap("bin_path", Config)['binbase'] + ConfigSectionMap("bwa", Config)['bwa_bin'],file, reverse_file, output_folder, reference, split_field, first_part, files_to_delete, logger, Config, type, command_list)
                out_sam = files_to_delete[0]
                command_list, files_to_delete = samtobam(out_sam, output_folder, analysis, files_to_delete, logger, Config, command_list, files_to_delete)
                out_bam = files_to_delete[1]
                command_list, files_to_delete = sort_bam(out_bam, output_folder, analysis, logger, Config, command_list, files_to_delete)
                out_sort_bam = files_to_delete[2]
                command_list = index_bam(out_sort_bam, output_folder, logger, Config, command_list, files_to_delete)
                command_list, gatk_depth_of_coverage_file = gatk_DepthOfCoverage(out_sorted_bam, output_folder, analysis, reference, logger, Config, command_list)
                command_list = flagstat(out_sorted_bam, coverage_depth_directory, analysis, logger, Config, command_list)
                coverage_depth_cmd = ""

                for i in command_list:
                    coverage_depth_cmd = coverage_depth_cmd + i + "\n"
                keep_logging("The coverage Depth commands for file %s are:\n", "The coverage Depth commands for file %s are:\n", logger, 'info')
                keep_logging(coverage_depth_cmd, coverage_depth_cmd, logger, 'debug')


                if cluster == "cluster":
                    generate_cluster_jobs(coverage_depth_cmd, file_prefix, Config, logger)
                else:
                    f3=open(file_prefix + '_commands.sh', 'w+')
                    f3.write(coverage_depth_cmd)

def prepare_readgroup(forward_read, logger):
    keep_logging('Preparing ReadGroup Info', 'Preparing ReadGroup Info', logger, 'info')
    samplename = os.path.basename(forward_read)
    if forward_read.endswith(".gz"):
        ###
        output = gzip.open(forward_read, 'rb')
        firstLine = output.readline()
        split_field = re.split(r":",firstLine)
        id_name = split_field[1]
        id_name = id_name.strip()
        split_field = "\"" + "@RG" + "\\tID:" + split_field[1] + "\\tSM:" + samplename + "\\tLB:1\\tPL:Illumina" + "\""
        return split_field

    elif forward_read.endswith(".fastq"):
        ###
        output = open(forward_read, 'r')
        firstLine = output.readline()
        split_field = re.split(r":",firstLine)
        split_field = "\"" + "@RG" + "\\tID:" + split_field[1] + "\\tSM:" + samplename + "\\tLB:1\\tPL:Illumina" + "\""
        return split_field


    elif forward_read.endswith(".fq"):
        ###
        output = open(forward_read, 'r')
        firstLine = output.readline()
        split_field = re.split(r":",firstLine)
        split_field = "\"" + "@RG" + "\\tID:" + split_field[1] + "\\tSM:" + samplename + "\\tLB:1\\tPL:Illumina" + "\""
        return split_field
