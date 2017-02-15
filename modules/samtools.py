__author__ = 'alipirani'
import os
from modules.log_modules import keep_logging
from modules.logging_subprocess import *
from config_settings import ConfigSectionMap

############################################################### SAM to BAM conversion ####################################################################################################
def samtobam(out_sam, out_path, analysis, files_to_delete, logger, Config, command_list):
    base_cmd = ConfigSectionMap("bin_path", Config)['binbase'] + "/" + ConfigSectionMap("samtools", Config)['samtools_bin'] + "/" + ConfigSectionMap("samtools", Config)['base_cmd']
    cmd = "%s view -Sb %s > %s/%s_aln.bam" % (base_cmd, out_sam, out_path, analysis)
    keep_logging('SAM to BAM Conversion', 'SAM to BAM Conversion', logger, 'info')
    keep_logging(cmd, cmd, logger, 'debug')
    try:
        command_list.append(cmd)
    except sp.CalledProcessError:
        keep_logging('Error in appending command to command list at SAM-to-BAM Conversion step. Exiting.', 'Error in appending command to command list at SAM-to-BAM Conversion step. Exiting.', logger, 'exception')
        sys.exit(1)
    out_bam = "%s/%s_aln.bam" % (out_path, analysis)
    files_to_delete.append(out_bam)
    return command_list, files_to_delete

############################################################### END: SAM to BAM conversion ###############################################################################################

############################################################### BAM Sorting ##############################################################################################################
def sort_bam(out_bam, out_path, analysis, logger, Config, command_list, files_to_delete):
    base_cmd = ConfigSectionMap("bin_path", Config)['binbase'] + "/" + ConfigSectionMap("samtools", Config)['samtools_bin'] + "/" + ConfigSectionMap("samtools", Config)['base_cmd']
    cmd = "%s sort %s %s/%s_aln_sort" % (base_cmd, out_bam, out_path, analysis)
    keep_logging('Sorting BAM file', 'Sorting BAM file', logger, 'info')
    keep_logging(cmd, cmd, logger, 'debug')
    try:
        command_list.append(cmd)
    except sp.CalledProcessError:
        keep_logging('Error in appending command to command list at BAM Sorting step. Exiting.', 'Error in appending command to command list at BAM Sorting step. Exiting.', logger, 'exception')
        sys.exit(1)
    sort_bam = "%s/%s_aln_sort.bam" % (out_path, analysis)
    files_to_delete.append(sort_bam)
    return command_list, files_to_delete
############################################################### END: BAM Sorting #########################################################################################################

############################################################### BAM Indexing ##############################################################################################################
def index_bam(out_sort_bam, out_path, logger, Config, command_list, files_to_delete):
    base_cmd = ConfigSectionMap("bin_path", Config)['binbase'] + "/" + ConfigSectionMap("samtools", Config)['samtools_bin'] + "/" + ConfigSectionMap("samtools", Config)['base_cmd']
    cmd = "%s index %s" % (base_cmd, out_sort_bam)
    keep_logging(cmd, cmd, logger, 'info')
    try:
        command_list.append(cmd)
    except sp.CalledProcessError:
        keep_logging('Error in appending command to command list at Samtools Indexing step. Exiting.', 'Error in appending command to command list at Samtools Indexing step. Exiting.', logger, 'exception')
        sys.exit(1)
    return command_list
############################################################### END: BAM Sorting ##########################################################################################################

# ############################################################### Reference FAI Indexing ####################################################################################################
# def ref_fai_index(reference):
#     cmd = "%s faidx %s" % (base_cmd, reference)
#     print "\nRunning:\n [%s] \n" % cmd
#     os.system(cmd)
# ############################################################### END: Reference FAI Indexing ###############################################################################################
#
#
############################################################################## Alignment Statistics: Flagstat ###############################################################################
def flagstat(out_sorted_bam, out_path, analysis, logger, Config, command_list):
    base_cmd = ConfigSectionMap("bin_path", Config)['binbase'] + "/" + ConfigSectionMap("samtools", Config)['samtools_bin'] + "/" + ConfigSectionMap("samtools", Config)['base_cmd']
    cmd = "%s flagstat %s > %s/%s_alignment_stats" % (base_cmd, out_sorted_bam, out_path, analysis)
    keep_logging(cmd, cmd, logger, 'debug')
    try:
        command_list.append(cmd)
    except sp.CalledProcessError:
        keep_logging('Error in appending command to command list at Samtools Alignment Stats step. Exiting.', 'Error in appending command to command list at Samtools Alignment Stats step. Exiting.', logger, 'exception')
        sys.exit(1)
    return command_list
############################################################################## END: Alignment Statistics #####################################################################################




