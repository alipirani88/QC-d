__author__ = 'alipirani'
import os
from modules.log_modules import keep_logging
from logging_subprocess import *




def align_bwa(base_cmd,forward_clean, reverse_clean, out_path, reference, split_field, analysis, files_to_delete, logger, Config, type, command_list):
    if type == "PE":
        cmd = "%s mem -M -R %s -t 8 %s %s %s > %s/%s_aln.sam" % (base_cmd,split_field, reference, forward_clean, reverse_clean, out_path, analysis)
    else:
        cmd = "%s mem -M -R %s -t 8 %s %s > %s/%s_aln.sam" % (base_cmd,split_field, reference, forward_clean, out_path, analysis)
    keep_logging(cmd, cmd, logger, 'debug')
    try:
        #call(cmd, logger)
        command_list.append(cmd)
    except sp.CalledProcessError:
        keep_logging('Error in Appending command to command list. Exiting.', 'Error in Appending command to command list. Exiting.', logger, 'exception')
        sys.exit(1)
    out_sam = "%s/%s_aln.sam" % (out_path, analysis)
    files_to_delete.append(out_sam)
    return command_list, files_to_delete