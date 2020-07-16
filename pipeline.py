__author__ = 'alipirani'

"""Declaring required python modules"""
import argparse
import ConfigParser
import subprocess
import re
import os
import sys
import errno
import glob
from modules.log_modules import  *
from logging_subprocess import *
from config_settings import ConfigSectionMap
from modules.coverage import *
from modules.quality import *
from modules.multiqc import *
from modules.screen_contamination import *
from modules.kraken_contamination import *
from datetime import datetime
from argparse import RawTextHelpFormatter
from modules.coverage_depth import *
from modules.kraken_report import *
from modules.mlst import mlst
from modules.summary import *

""" Command line argument parsing """
def parser():
    parser = argparse.ArgumentParser(description='QC\'d - Quality control and Contamination Detection Workflow', formatter_class=RawTextHelpFormatter)
    required = parser.add_argument_group('Required arguments')
    optional = parser.add_argument_group('Optional arguments')
    required.add_argument('-samples', action='store', dest="samples", help='A file containing filenames of forward-paired end or single-end reads. One Sample per line')
    optional.add_argument('-config', action='store', dest="config", help='Path to Config file, Make sure to check config settings before running pipeline.\nNote: Set Kraken database path under [kraken] config section')
    required.add_argument('-dir', action='store', dest="directory", help='Path to Sequencing Reads Data directory. NOTE: Provide full/absolute path.')
    required.add_argument('-analysis', action='store', dest="analysis_names", help='comma-seperated analysis names to run\n[Options: coverage,quality,kraken_contamination,kraken_report,coverage_depth].\nExample: "-analysis coverage,quality" - This will estimate coverage and quality for the samples given in -samples')
    required.add_argument('-o', action='store', dest="output_folder", help='Output Folder Path ending with output directory name to save the results.\nCreates a new output directory path if it doesn\'t exist.')
    required.add_argument('-type', action='store', dest='type', help='Type of analysis: SE or PE')
    optional.add_argument('-cluster', action='store', dest='cluster', help='Run in one of the two modes. \nDefault is local for coverage and fastqc analysis.\nFor all the other analysis Default is cluster. The pipeline prefers cluster mode to generate SLURM/PBS jobs.\nModes: local or cluster')
    required.add_argument('-genome_size', action='store', dest='size', help='Estimated Genome Size - to be used for estimating coverage in downsampling steps')
    required.add_argument('-prefix', action='store', dest='prefix', help='Use this prefix to save the results files')
    optional.add_argument('-reference', action='store', dest='reference', help='Reference genome to use for calculating GATK Depth of coverage.')
    optional.add_argument('-downsample', action='store', dest="downsample",
                          help='yes/no: Downsample Reads data to default depth of 100X')
    optional.add_argument('-scheduler', action='store', dest="scheduler",
                          help='\nType of Scheduler for generating Kraken cluster jobs: PBS, SLURM, LOCAL')
    optional.add_argument('-dryrun', action='store_true', dest="dryrun",
                          help='Perform a trial run without submitting cluster jobs',
                          required=False)
    optional.add_argument('-mlst_db', action='store', dest="mlst_db",
                          help='Ariba MLST database path',
                          required=False)
    return parser

""" Main Pipeline """
def pipeline(args, logger, Config, output_folder, prefix, reference):
    keep_logging('\nSTART: Pipeline', 'START: Pipeline', logger, 'info')

    """ Check Subroutines and create logger object: Arguments, Input files, Reference Index"""
    keep_logging('Checking Dependencies...', 'Checking Dependencies', logger, 'info')

    """ Check java availability """
    java_check()

    """ Check if the input file exists """
    with open(args.samples) as fp:
        for line in fp:
            line = line.strip()
            line = args.directory + "/" + line
            filenames_array.append(line)
            if args.type != "PE":
                reverse_raw = "None"
                file_exists(line, reverse_raw)
            else:
                #reverse_raw = args.directory + "/" + reverse_raw
                file_exists(line, line)

    keep_logging('Total no. of Samples %s' % len(filenames_array), 'Total no. of Samples %s' % len(filenames_array), logger, 'info')

    """ Start the pipeline: """
    analysis_list = args.analysis_names.split(',')

    keep_logging('Running analysis - %s' % args.analysis_names, 'Running analysis - %s' % args.analysis_names,
                 logger, 'info')

    """ Copy filenames to output folder """
    cp_cmd = "cp %s %s" % (args.samples, output_folder)
    os.system(cp_cmd)

    """ Set Default cluster mode"""
    if args.cluster:
        cluster = args.cluster
    else:
        cluster = "local"

    """ Start Specific analysis based on analysis list """
    for analysis in analysis_list:
        if analysis == "coverage":
            keep_logging("Step: Calculating Coverage...\n", "Calculating Coverage", logger, 'info')
            coverage(filenames_array, Config, logger, output_folder, args.type, args.samples, args.size, prefix)
        elif analysis == "quality":
            keep_logging("Step: Analysing Fastqc Quality...\n", "Analysing Fastqc Quality...", logger, 'info')
            fastqc_main_directory = args.output_folder + "/%s_Fastqc" % args.prefix
            make_sure_path_exists(fastqc_main_directory)
            fastqc_forward_directory = fastqc_main_directory + "/%s_Forward" % args.prefix
            make_sure_path_exists(fastqc_forward_directory)
            fastqc_reverse_directory = fastqc_main_directory + "/%s_Reverse" % args.prefix
            make_sure_path_exists(fastqc_reverse_directory)
            Multiqc_reports_directory = args.output_folder + "/%s_Multiqc_reports" % args.prefix
            make_sure_path_exists(Multiqc_reports_directory)
            quality(filenames_array, Config, logger, output_folder, args.type, args.samples, fastqc_forward_directory, fastqc_reverse_directory)
            multiqc(fastqc_forward_directory, "%s_Forward_fastqc" % args.prefix, Config, logger, Multiqc_reports_directory)
            multiqc(fastqc_reverse_directory, "%s_Reverse_fastqc" % args.prefix, Config, logger, Multiqc_reports_directory)
        elif analysis == "screen_contamination":
            keep_logging("Step: Screening Fastq reads against Reference Database...\n", "Screening Fastq reads against Reference Database...", logger, 'info')
            fastq_screen_directory = args.output_folder + "/%s_Fastqc_screen" % args.prefix
            make_sure_path_exists(fastq_screen_directory)
            screen_contamination(filenames_array, Config, logger, output_folder, args.type, args.samples, fastq_screen_directory, cluster)
            Multiqc_reports_directory = args.output_folder + "/%s_Multiqc_reports" % args.prefix
            make_sure_path_exists(Multiqc_reports_directory)
            multiqc(fastq_screen_directory, "%s_Fastq_screen" % args.prefix, Config, logger, Multiqc_reports_directory)
            keep_logging('MultiQC Report of FastQC results can be found in - %s\n' % Multiqc_reports_directory, 'MultiQC Report of FastQC results can be found in - %s\n' % Multiqc_reports_directory, logger, 'info')
        elif analysis == "kraken_contamination":
            keep_logging("Step: Running Kraken on Input reads...\n", "Running Kraken on Input reads...", logger, 'info')
            kraken_directory = args.output_folder + "/%s_Kraken_results" % args.prefix
            make_sure_path_exists(kraken_directory)
            kraken_contamination(filenames_array, Config, logger, output_folder, args.type, args.samples, kraken_directory, cluster, args.downsample, args.scheduler, args.size, args.dryrun)
        elif analysis == "kraken_report":
            keep_logging("Step: Generating Kraken report on Kraken Results...\n", "Generating Kraken report on Kraken Results...", logger, 'info')
            kraken_directory = args.output_folder + "/%s_Kraken_results" % args.prefix
            make_sure_path_exists(kraken_directory)
            kraken_report(filenames_array, Config, logger, output_folder, args.type, args.samples, kraken_directory, cluster, args.scheduler)
        elif analysis == "coverage_depth":
            keep_logging("Step: Running Coverage Depth analysis on Input reads...\n", "Running Coverage Depth analysis on Input reads...", logger, 'info')
            coverage_depth_directory = args.output_folder + "/%s_Coverage_depth" % args.prefix
            make_sure_path_exists(coverage_depth_directory)
            coverage_depth_analysis(filenames_array, Config, logger, output_folder, args.type, args.samples, coverage_depth_directory, cluster, reference, args.scheduler)
        elif analysis == "mlst":
            keep_logging("Step: Running Ariba MLST sequence typing on Input reads...\n", "Running MLST sequence typing on Input reads...", logger, 'info')
            if args.mlst_db:
                mlstdb = args.mlst_db
            else:
                mlstdb = ConfigSectionMap("ariba", Config)['mlst_db_path']
            keep_logging(
                '',
                "Using Ariba MLST Database from this path - %s" % mlstdb,
                logger, 'debug')
            mlst_directory = args.output_folder + "/%s_MLST_results" % args.prefix
            make_sure_path_exists(mlst_directory)
            mlst(filenames_array, Config, logger, mlst_directory, args.type, args.samples, mlst_directory, cluster, args.scheduler, mlstdb)
        elif analysis == "summary":
            keep_logging('', "Generating Summary report for QC'd analysis - %s" % args.prefix, logger, 'debug')
            summary(filenames_array, Config, logger, args.prefix, output_folder)
            keep_logging("Summary report - %s/%s_summary.tsv" % (output_folder, prefix), "Summary report - %s/%s_summary.tsv" % (output_folder, prefix), logger, 'info')


""" Check Subroutines """

""" Usage Message: """
def usage():
    print "Usage: pipeline.py [-h] [-samples FILENAMES] [-dir DIRECTORY of Fastq files] [-config CONFIG] [-analysis ANALYSIS_NAME1,ANALYSIS_NAME2,...] [-o OUTPUT_FOLDER] [-type TYPE]\n"

""" Check Java Availability: """
def java_check():
    print "Checking Java Availability....\n"
    jd = subprocess.check_output(["java", "-version"], stderr=subprocess.STDOUT)
    if len(jd) < 1:
        print "Unable to find a java runtime environment. The pipeline requires java 6 or later.\n"
    else:
        print "Java Availability Check completed ...\n" + jd

""" Make sure input raw reads files exists at given location. """
def file_exists(path1, path2):
    if path1:
        if not os.path.isfile(path1):
            file_basename = os.path.basename(path1)
            print "The input file " + file_basename + " does not exists. \nPlease provide another file.\n"
            exit()
    if path2 != "None":
        if not os.path.isfile(path2):
            file_basename = os.path.basename(path2)
            print "The input file " + file_basename + " does not exists. \nPlease provide another file.\n"
            exit()

""" Make sure the output folder exists or create at given path """
def make_sure_path_exists(out_path):
    try:
        os.makedirs(out_path)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            print "Errors in output folder path! please change the output path or analysis name\n"
            exit()

""" Start of Main Method/Pipeline """
if __name__ == '__main__':
    start_time = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    start_time_2 = datetime.now()
    log_unique_time = datetime.now().strftime('%Y_%m_%d_%H_%M_%S')

    args = parser().parse_args()

    global config_file
    global filenames_array
    global Config
    global logger
    analysis_string = args.analysis_names.replace(',', '_')

    filenames_array = []

    # Set config file paths
    if args.config:
        config_file = args.config
    else:
        config_file = os.path.dirname(os.path.abspath(__file__)) + "/config"
    Config = ConfigParser.ConfigParser()
    Config.read(config_file)

    logger = generate_logger(args.output_folder, analysis_string, log_unique_time)




    # Set output directory paths
    if args.output_folder != '':
        args.output_folder += '/'
        make_sure_path_exists(args.output_folder)
    if args.output_folder != '':
        args.output_folder += '/'

    # Set reference genome path to map the samples and calculate coverage depth
    if "coverage_depth" in args.analysis_names:
        try:
            reference = ConfigSectionMap(args.reference, Config)['ref_path'] + "/" + ConfigSectionMap(args.reference, Config)['ref_name']
        except OSError as exception:
            if exception.errno != errno.EEXIST:
                print "Please provide reference genome name or Check the reference genome path in config file.\n"
                exit()
    else:
        reference = "NONE"

    # Main Workflow
    pipeline(args, logger, Config, args.output_folder, args.prefix, reference)

    keep_logging('End: Pipeline\n', 'End: Pipeline', logger, 'info')
    time_taken = datetime.now() - start_time_2
    keep_logging('Total Time taken: {}'.format(time_taken), 'Total Time taken: {}'.format(time_taken), logger, 'info')