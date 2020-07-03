__author__ = 'alipirani'
import os
from modules.log_modules import keep_logging
from logging_subprocess import *
from config_settings import ConfigSectionMap


def gatk_DepthOfCoverage(out_sorted_bam, out_path, analysis_name, reference, logger, Config, command_list):
    base_cmd = os.path.dirname(os.path.dirname(os.path.abspath(__file__))) + "/" + ConfigSectionMap("gatk", Config)['base_cmd']
    cmd = "java -jar %s -T DepthOfCoverage -R %s -o %s/%s_depth_of_coverage -I %s --summaryCoverageThreshold 1 --summaryCoverageThreshold 5 --summaryCoverageThreshold 9 --summaryCoverageThreshold 10 --summaryCoverageThreshold 15 --summaryCoverageThreshold 20 --summaryCoverageThreshold 25 --minBaseQuality 15" % (base_cmd, reference, out_path, analysis_name, out_sorted_bam)
    keep_logging('', cmd, logger, 'debug')
    try:
        command_list.append(cmd)
    except sp.CalledProcessError:
        keep_logging('Error in GATK Depth of Coverage step. Exiting.', 'Error in GATK Depth of Coverage step. Exiting.', logger, 'exception')
        sys.exit(1)
    gatk_depth_of_coverage_file = "%s/%s_depth_of_coverage" % (out_path, analysis_name)
    keep_logging('', 'GATK Depth of Coverage file: {}'.format(gatk_depth_of_coverage_file), logger, 'debug')
    return command_list, gatk_depth_of_coverage_file