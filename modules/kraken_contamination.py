__author__ = 'alipirani'
import os
import subprocess
from modules.log_modules import keep_logging
from modules.logging_subprocess import *
from config_settings import ConfigSectionMap
from modules.krona_visualization import *
from modules.generate_cluster_jobs import *
from modules.run_parallel import *

def kraken_contamination(filenames_array, Config, logger, output_folder, type, samples, kraken_directory, cluster):
    parallel_local_cmds = []
    parallel_local_cmds_krona = []
    cmd = ""
    if type == "PE":
        for file in filenames_array:
            file_prefix = kraken_directory + "/" + os.path.basename(file)[0:20]
            if file.endswith('.gz'):
                #kraken_cmd = "%s/%s/%s --fastq-input --gzip-compressed --unclassified-out %s_unclassified.txt --db %s --output %s --paired %s %s" % (ConfigSectionMap("bin_path", Config)['binbase'], ConfigSectionMap("kraken", Config)['kraken_bin'], ConfigSectionMap("kraken", Config)['base_cmd'], file_prefix, ConfigSectionMap("kraken", Config)['db_path'], file_prefix, file, file.replace('_R1_', '_R2_'))
                kraken_cmd = "%s/%s/%s --fastq-input --gzip-compressed --unclassified-out %s_unclassified.txt --db %s --output %s_kraken %s" % (ConfigSectionMap("bin_path", Config)['binbase'], ConfigSectionMap("kraken", Config)['kraken_bin'], ConfigSectionMap("kraken", Config)['base_cmd'], file_prefix, ConfigSectionMap("kraken", Config)['db_path'], file_prefix, file)
                keep_logging(kraken_cmd, kraken_cmd, logger, 'debug')
                if cluster == "cluster":
                    cmd = cmd + "\n" + kraken_cmd
                    krona_cmd = krona_visualization(file_prefix, Config, logger, kraken_directory, cluster)
                    cmd = cmd + "\n" + krona_cmd
                    generate_cluster_jobs(cmd, file_prefix, Config, logger)
                elif cluster == "local":
                    call(kraken_cmd, logger)
                    krona_cmd = krona_visualization(file_prefix, Config, logger, kraken_directory, cluster)
                elif cluster == "parallel-local":
                    #parallel_local_cmds.append(kraken_cmd)
                    krona_cmd = krona_visualization(file_prefix, Config, logger, kraken_directory, cluster)
                    parallel_local_cmds_krona.append(krona_cmd)
            else:
                #kraken_cmd = "%s/%s/%s --fastq-input --gzip-compressed --unclassified-out %s_unclassified.txt --db %s --output %s --paired %s %s" % (ConfigSectionMap("bin_path", Config)['binbase'], ConfigSectionMap("kraken", Config)['kraken_bin'], ConfigSectionMap("kraken", Config)['base_cmd'], file_prefix, ConfigSectionMap("kraken", Config)['db_path'], file_prefix, file, file.replace('_R1_', '_R2_'))
                kraken_cmd = "%s/%s/%s --fastq-input --unclassified-out %s_unclassified.txt --db %s --output %s_kraken %s" % (ConfigSectionMap("bin_path", Config)['binbase'], ConfigSectionMap("kraken", Config)['kraken_bin'], ConfigSectionMap("kraken", Config)['base_cmd'], file_prefix, ConfigSectionMap("kraken", Config)['db_path'], file_prefix, file)
                keep_logging(kraken_cmd, kraken_cmd, logger, 'debug')
                if cluster == "cluster":
                    cmd = cmd + "\n" + kraken_cmd
                    krona_cmd = krona_visualization(file_prefix, Config, logger, kraken_directory, cluster)
                    cmd = cmd + "\n" + krona_cmd
                    generate_cluster_jobs(cmd, file_prefix, Config, logger)
                elif cluster == "local":
                    call(kraken_cmd, logger)
                    krona_cmd = krona_visualization(file_prefix, Config, logger, kraken_directory, cluster)
                elif cluster == "parallel-local":
                    parallel_local_cmds.append(kraken_cmd)
                    krona_cmd = krona_visualization(file_prefix, Config, logger, kraken_directory, cluster)
                    parallel_local_cmds.append(krona_cmd)
    elif type == "SE":
        for file in filenames_array:
            file_prefix = kraken_directory + "/" + os.path.basename(file)[0:20]
            if file.endswith('.gz'):
                #kraken_cmd = "%s/%s/%s --fastq-input --gzip-compressed --unclassified-out %s_unclassified.txt --db %s --output %s %s" % (ConfigSectionMap("bin_path", Config)['binbase'], ConfigSectionMap("kraken", Config)['kraken_bin'], ConfigSectionMap("kraken", Config)['base_cmd'], file_prefix, ConfigSectionMap("kraken", Config)['db_path'], file_prefix, file)
                kraken_cmd = "%s/%s/%s --fastq-input --gzip-compressed --unclassified-out %s_unclassified.txt --db %s --output %s_kraken %s" % (ConfigSectionMap("bin_path", Config)['binbase'], ConfigSectionMap("kraken", Config)['kraken_bin'], ConfigSectionMap("kraken", Config)['base_cmd'], file_prefix, ConfigSectionMap("kraken", Config)['db_path'], file_prefix, file)
                keep_logging(kraken_cmd, kraken_cmd, logger, 'debug')
                if cluster == "cluster":
                    cmd = cmd + "\n" + kraken_cmd
                    krona_cmd = krona_visualization(file_prefix, Config, logger, kraken_directory, cluster)
                    cmd = cmd + "\n" + krona_cmd
                    generate_cluster_jobs(cmd, file_prefix, Config, logger)
                elif cluster == "local":
                    call(kraken_cmd, logger)
                    krona_cmd = krona_visualization(file_prefix, Config, logger, kraken_directory, cluster)
                elif cluster == "parallel-local":
                    parallel_local_cmds.append(kraken_cmd)
                    krona_cmd = krona_visualization(file_prefix, Config, logger, kraken_directory, cluster)
                    parallel_local_cmds.append(krona_cmd)
            else:
                #kraken_cmd = "%s/%s/%s --fastq-input --gzip-compressed --unclassified-out %s_unclassified.txt --db %s --output %s %s" % (ConfigSectionMap("bin_path", Config)['binbase'], ConfigSectionMap("kraken", Config)['kraken_bin'], ConfigSectionMap("kraken", Config)['base_cmd'], file_prefix, ConfigSectionMap("kraken", Config)['db_path'], file_prefix, file)
                kraken_cmd = "%s/%s/%s --fastq-input --unclassified-out %s_unclassified.txt --db %s --output %s_kraken %s" % (ConfigSectionMap("bin_path", Config)['binbase'], ConfigSectionMap("kraken", Config)['kraken_bin'], ConfigSectionMap("kraken", Config)['base_cmd'], file_prefix, ConfigSectionMap("kraken", Config)['db_path'], file_prefix, file)
                keep_logging(kraken_cmd, kraken_cmd, logger, 'debug')
                if cluster == "cluster":
                    cmd = cmd + "\n" + kraken_cmd
                    krona_cmd = krona_visualization(file_prefix, Config, logger, kraken_directory, cluster)
                    cmd = cmd + "\n" + krona_cmd
                    generate_cluster_jobs(cmd, file_prefix, Config, logger)
                elif cluster == "local":
                    call(kraken_cmd, logger)
                    krona_cmd = krona_visualization(file_prefix, Config, logger, kraken_directory, cluster)
                elif cluster == "parallel-local":
                    parallel_local_cmds.append(kraken_cmd)
                    krona_cmd = krona_visualization(file_prefix, Config, logger, kraken_directory, cluster)
                    parallel_local_cmds.append(krona_cmd)
    if cluster == "parallel-local":
        complete = run_parallel(parallel_local_cmds)
        complete = run_parallel(parallel_local_cmds_krona)