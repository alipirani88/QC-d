__author__ = 'alipirani'
import os
import subprocess
from modules.log_modules import keep_logging
from logging_subprocess import *
from config_settings import ConfigSectionMap
from modules.krona_visualization import *
from modules.generate_cluster_jobs import *
from modules.run_parallel import *


def downsample_reads(file, file2, logger):


    # Run Mash to estimate Genome size
    keep_logging('Running: /nfs/esnitkin/bin_group/variant_calling_bin/mash sketch -o /tmp/sketch_out -k 32 -m 3 -r %s >& /tmp/sketch_stdout' % file,
                 'Running: /nfs/esnitkin/bin_group/variant_calling_bin/mash sketch -o /tmp/sketch_out -k 32 -m 3 -r %s >& /tmp/sketch_stdout' % file, logger, 'info')

    mash_cmd = "/nfs/esnitkin/bin_group/variant_calling_bin/mash sketch -o /tmp/sketch_out -k 32 -m 3 -r %s >& /tmp/sketch_stdout" % file

    keep_logging('Running: %s' % mash_cmd,
                 'Running: %s' % mash_cmd, logger, 'info')


    try:
        call(mash_cmd, logger)
    except sp.CalledProcessError:
        keep_logging('Error running Mash for estimating genome size.', 'Error running Mash for estimating genome size', logger, 'exception')
        sys.exit(1)

    with open("/tmp/sketch_stdout", 'rU') as file_open:
        for line in file_open:
            if line.startswith('Estimated genome size:'):
                gsize = float(line.split(': ')[1].strip())
            if line.startswith('Estimated coverage:'):
                est_cov = float(line.split(': ')[1].strip())
    file_open.close()


    keep_logging('Estimated Genome Size from Mash Sketch: %s' % gsize,
                 'Estimated Genome Size from Mash Sketch: %s' % gsize, logger, 'info')

    # Extract basic fastq reads stats with seqtk
    seqtk_check = "/nfs/esnitkin/bin_group/seqtk/seqtk fqchk -q3 %s > /tmp/%s_fastqchk.txt" % (file, os.path.basename(file))

    keep_logging('Running seqtk to extract Fastq statistics: %s' % seqtk_check,
                 'Running seqtk to extract Fastq statistics: %s' % seqtk_check, logger, 'info')


    try:
        call(seqtk_check, logger)
    except sp.CalledProcessError:
        keep_logging('Error running seqtk for extracting fastq statistics.', 'Error running seqtk for extracting fastq statistics.', logger, 'exception')
        sys.exit(1)


    with open("/tmp/%s_fastqchk.txt" % os.path.basename(file), 'rU') as file_open:
        for line in file_open:
            if line.startswith('min_len'):
                line_split = line.split(';')
                min_len = line_split[0].split(': ')[1]
                max_len = line_split[1].split(': ')[1]
                avg_len = line_split[2].split(': ')[1]
            if line.startswith('ALL'):
                line_split = line.split('\t')
                total_bases = int(line_split[1]) * 2
    file_open.close()


    keep_logging('Average Read Length: %s' % avg_len,
                 'Average Read Length: %s' % avg_len, logger, 'info')

    keep_logging('Total number of bases in fastq: %s' % total_bases,
                 'Total number of bases in fastq: %s' % total_bases, logger, 'info')

    # Calculate original depth and check if it needs to be downsampled to a default coverage.
    ori_coverage_depth = int(total_bases / gsize)

    keep_logging('Original Covarage Depth: %s x' % ori_coverage_depth,
                 'Original Covarage Depth: %s x' % ori_coverage_depth, logger, 'info')

    proc = sp.Popen(["nproc"], stdout=sp.PIPE, shell=True)
    (nproc, err) = proc.communicate()
    nproc = nproc.strip()

    if ori_coverage_depth > 100:
        # Downsample to 100
        factor = float(100 / float(ori_coverage_depth))
        print factor
        r1_sub = "/tmp/%s" % os.path.basename(file)

        # Downsample using seqtk
        try:
            keep_logging("/nfs/esnitkin/bin_group/seqtk/seqtk sample %s %s | pigz --fast -c -p %s > /tmp/%s" % (
                file, factor, nproc, os.path.basename(file)),
                         "/nfs/esnitkin/bin_group/seqtk/seqtk sample %s %s | pigz --fast -c -p %s > /tmp/%s" % (
                file, factor, nproc, os.path.basename(file)), logger, 'info')
            # call("/nfs/esnitkin/bin_group/seqtk/seqtk sample %s %s | pigz --fast -c -p %s > /tmp/%s" % (
            #     file, factor, nproc, os.path.basename(file)), logger)
            seqtk_downsample = "/nfs/esnitkin/bin_group/seqtk/seqtk sample %s %s | pigz --fast -c -p %s > /tmp/%s" % (
                file, factor, nproc, os.path.basename(file))
        except sp.CalledProcessError:
            keep_logging('Error running seqtk for downsampling raw fastq reads.',
                         'Error running seqtk for downsampling raw fastq reads.', logger, 'info')
            sys.exit(1)

        if file2:
            r2_sub = "/tmp/%s" % os.path.basename(file2)

            try:
                keep_logging("/nfs/esnitkin/bin_group/seqtk/seqtk sample %s %s | pigz --fast -c -p %s > /tmp/%s" % (
                    file2, factor, nproc, os.path.basename(file2)),
                             "/nfs/esnitkin/bin_group/seqtk/seqtk sample %s %s | pigz --fast -c -p %s > /tmp/%s" % (
                    file2, factor, nproc, os.path.basename(file2)), logger, 'info')
                # call("/nfs/esnitkin/bin_group/seqtk/seqtk sample %s %s | pigz --fast -c -p %s > /tmp/%s" % (
                #     file2, factor, nproc, os.path.basename(file2)), logger)
            except sp.CalledProcessError:
                keep_logging('Error running seqtk for downsampling raw fastq reads.',
                             'Error running seqtk for downsampling raw fastq reads.', logger, 'exception')
                sys.exit(1)
        else:
            r2_sub = "None"

    elif ori_coverage_depth < 100:
        r1_sub = file
        r2_sub = file2
        seqtk_downsample = ""

    return r1_sub, r2_sub, seqtk_downsample


def kraken_contamination(filenames_array, Config, logger, output_folder, type, samples, kraken_directory, cluster, downsample, scheduler):
    parallel_local_cmds = []
    parallel_local_cmds_krona = []
    cmd = ""
    for file in filenames_array:
        file_prefix = kraken_directory + "/" + os.path.basename(file)[0:20]
        file2 = file.replace('_R1_', '_R2_')
        if downsample == "yes":
            read1, read2, seqtk_downsample = downsample_reads(file, file2, logger)
            file = read1
            keep_logging("Using downsampled reads - %s" % file, "Using downsampled reads - %s" % file, logger, 'info')

        kraken_cmd = seqtk_downsample
        if file.endswith('.gz'):
            file = "/tmp/%s" % os.path.basename(file)
            kraken_cmd = kraken_cmd + '\n' + "%s/%s/%s --quick --fastq-input --gzip-compressed --unclassified-out %s_unclassified.txt --db %s --output %s_kraken %s" % (ConfigSectionMap("bin_path", Config)['binbase'], ConfigSectionMap("kraken", Config)['kraken_bin'], ConfigSectionMap("kraken", Config)['base_cmd'], file_prefix, ConfigSectionMap("kraken", Config)['db_path'], file_prefix, file)
            keep_logging(kraken_cmd, kraken_cmd, logger, 'info')
            # kraken_cmd = "%s/%s/%s --quick --gzip-compressed --unclassified-out %s_unclassified.txt --db %s --output %s_kraken --report %s_kraken_report.txt %s" % (
            # ConfigSectionMap("bin_path", Config)['binbase'], ConfigSectionMap("kraken2", Config)['kraken_bin'],
            # ConfigSectionMap("kraken2", Config)['base_cmd'], file_prefix, ConfigSectionMap("kraken2", Config)['db_path'],
            # file_prefix, file_prefix, file)
            keep_logging(kraken_cmd, kraken_cmd, logger, 'info')
            kraken_out = file_prefix + "_kraken"
            report_cmd = "%s/%s/kraken-report --db %s %s > %s_report.txt" % (
            ConfigSectionMap("bin_path", Config)['binbase'], ConfigSectionMap("kraken", Config)['kraken_bin'],
            ConfigSectionMap("kraken", Config)['db_path'], kraken_out, kraken_out)
            kraken_commands =  kraken_cmd + "\n" + report_cmd
            if cluster == "cluster":
                cmd = cmd + "\n" + kraken_cmd
                krona_cmd = krona_visualization(file_prefix, Config, logger, kraken_directory, cluster)
                cmd = cmd + "\n" + krona_cmd
                generate_cluster_jobs(kraken_commands, file_prefix, scheduler, Config, logger)
            elif cluster == "local":
                call(kraken_cmd, logger)
                #keep_logging(kraken_cmd, kraken_cmd, logger, 'info')
                krona_cmd = krona_visualization(file_prefix, Config, logger, kraken_directory, cluster)
                keep_logging(krona_cmd, krona_cmd, logger, 'info')
                parallel_local_cmds.append(kraken_cmd)
                parallel_local_cmds_krona.append(krona_cmd)
                call(krona_cmd, logger)
            elif cluster == "parallel-local":
                parallel_local_cmds.append(kraken_cmd)
                krona_cmd = krona_visualization(file_prefix, Config, logger, kraken_directory, cluster)
                parallel_local_cmds_krona.append(krona_cmd)
        else:
            kraken_cmd = "%s/%s/%s --quick --unclassified-out %s_unclassified.txt --db %s --output %s_kraken --report %s_kraken_report.txt %s" % (
                ConfigSectionMap("bin_path", Config)['binbase'], ConfigSectionMap("kraken2", Config)['kraken_bin'],
                ConfigSectionMap("kraken2", Config)['base_cmd'], file_prefix,
                ConfigSectionMap("kraken2", Config)['db_path'],
                file_prefix, file_prefix, file)
            if cluster == "cluster":
                cmd = cmd + "\n" + kraken_cmd
                krona_cmd = krona_visualization(file_prefix, Config, logger, kraken_directory, cluster)
                cmd = cmd + "\n" + krona_cmd
                generate_cluster_jobs(cmd, file_prefix, scheduler, Config, logger)
            elif cluster == "local":
                call(kraken_cmd, logger)
                #keep_logging(kraken_cmd, kraken_cmd, logger, 'info')
                krona_cmd = krona_visualization(file_prefix, Config, logger, kraken_directory, cluster)
                keep_logging(krona_cmd, krona_cmd, logger, 'info')
                parallel_local_cmds.append(kraken_cmd)
                parallel_local_cmds_krona.append(krona_cmd)
                call(krona_cmd, logger)
            elif cluster == "parallel-local":
                parallel_local_cmds.append(kraken_cmd)
                krona_cmd = krona_visualization(file_prefix, Config, logger, kraken_directory, cluster)
                parallel_local_cmds_krona.append(krona_cmd)

    if cluster == "parallel-local":
        if len(parallel_local_cmds) > 1:
            complete = run_parallel(parallel_local_cmds)
        else:
            call(kraken_cmd, logger)
        if len(parallel_local_cmds_krona) > 1:
            complete = run_parallel(parallel_local_cmds_krona)
        else:
            call(krona_cmd, logger)
    elif cluster == "cluster":
        print "Generated jobs in the kraken results folder - %s" % kraken_directory
    elif cluster == "local":
        for i in parallel_local_cmds:
            print i



# def kraken_contamination(filenames_array, Config, logger, output_folder, type, samples, kraken_directory, cluster):
#     parallel_local_cmds = []
#     parallel_local_cmds_krona = []
#     cmd = ""
#     if type == "PE":
#         for file in filenames_array:
#             file_prefix = kraken_directory + "/" + os.path.basename(file)[0:20]
#             if file.endswith('.gz'):
#                 kraken_cmd = "%s/%s/%s --quick --fastq-input --gzip-compressed --unclassified-out %s_unclassified.txt --db %s --output %s_kraken %s" % (ConfigSectionMap("bin_path", Config)['binbase'], ConfigSectionMap("kraken", Config)['kraken_bin'], ConfigSectionMap("kraken", Config)['base_cmd'], file_prefix, ConfigSectionMap("kraken", Config)['db_path'], file_prefix, file)
#                 keep_logging('', kraken_cmd, logger, 'debug')
#                 if cluster == "cluster":
#                     cmd = cmd + "\n" + kraken_cmd
#                     krona_cmd = krona_visualization(file_prefix, Config, logger, kraken_directory, cluster)
#                     cmd = cmd + "\n" + krona_cmd
#                     generate_cluster_jobs(cmd, file_prefix, Config, logger)
#                 elif cluster == "local":
#                     call(kraken_cmd, logger)
#                     keep_logging(kraken_cmd, kraken_cmd, logger, 'info')
#                     #keep_logging(kraken_cmd, kraken_cmd, logger, 'info')
#                     krona_cmd = krona_visualization(file_prefix, Config, logger, kraken_directory, cluster)
#                     keep_logging(krona_cmd, krona_cmd, logger, 'info')
#                     parallel_local_cmds.append(kraken_cmd)
#                     parallel_local_cmds_krona.append(krona_cmd)
#                     call(krona_cmd, logger)
#                 elif cluster == "parallel-local":
#                     parallel_local_cmds.append(kraken_cmd)
#                     krona_cmd = krona_visualization(file_prefix, Config, logger, kraken_directory, cluster)
#                     parallel_local_cmds_krona.append(krona_cmd)
#             else:
#                 kraken_cmd = "%s/%s/%s --quick --fastq-input --unclassified-out %s_unclassified.txt --db %s --output %s_kraken %s" % (ConfigSectionMap("bin_path", Config)['binbase'], ConfigSectionMap("kraken", Config)['kraken_bin'], ConfigSectionMap("kraken", Config)['base_cmd'], file_prefix, ConfigSectionMap("kraken", Config)['db_path'], file_prefix, file)
#                 keep_logging('', kraken_cmd, logger, 'debug')
#                 if cluster == "cluster":
#                     cmd = cmd + "\n" + kraken_cmd
#                     krona_cmd = krona_visualization(file_prefix, Config, logger, kraken_directory, cluster)
#                     cmd = cmd + "\n" + krona_cmd
#                     generate_cluster_jobs(cmd, file_prefix, Config, logger)
#                 elif cluster == "local":
#                     #call(kraken_cmd, logger)
#                     krona_cmd = krona_visualization(file_prefix, Config, logger, kraken_directory, cluster)
#                 elif cluster == "parallel-local":
#                     parallel_local_cmds.append(kraken_cmd)
#                     krona_cmd = krona_visualization(file_prefix, Config, logger, kraken_directory, cluster)
#                     parallel_local_cmds_krona.append(krona_cmd)
#     elif type == "SE":
#         for file in filenames_array:
#             file_prefix = kraken_directory + "/" + os.path.basename(file)[0:20]
#             if file.endswith('.gz'):
#                 kraken_cmd = "%s/%s/%s --quick --fastq-input --gzip-compressed --unclassified-out %s_unclassified.txt --db %s --output %s_kraken %s" % (ConfigSectionMap("bin_path", Config)['binbase'], ConfigSectionMap("kraken", Config)['kraken_bin'], ConfigSectionMap("kraken", Config)['base_cmd'], file_prefix, ConfigSectionMap("kraken", Config)['db_path'], file_prefix, file)
#                 keep_logging('', kraken_cmd, logger, 'debug')
#                 if cluster == "cluster":
#                     cmd = cmd + "\n" + kraken_cmd
#                     krona_cmd = krona_visualization(file_prefix, Config, logger, kraken_directory, cluster)
#                     cmd = cmd + "\n" + krona_cmd
#                     generate_cluster_jobs(cmd, file_prefix, Config, logger)
#                 elif cluster == "local":
#                     call(kraken_cmd, logger)
#                     krona_cmd = krona_visualization(file_prefix, Config, logger, kraken_directory, cluster)
#                     call(krona_cmd)
#                 elif cluster == "parallel-local":
#                     parallel_local_cmds.append(kraken_cmd)
#                     krona_cmd = krona_visualization(file_prefix, Config, logger, kraken_directory, cluster)
#                     parallel_local_cmds_krona.append(krona_cmd)
#             else:
#                 kraken_cmd = "%s/%s/%s --fastq-input --unclassified-out %s_unclassified.txt --db %s --output %s_kraken %s" % (ConfigSectionMap("bin_path", Config)['binbase'], ConfigSectionMap("kraken", Config)['kraken_bin'], ConfigSectionMap("kraken", Config)['base_cmd'], file_prefix, ConfigSectionMap("kraken", Config)['db_path'], file_prefix, file)
#                 keep_logging('', kraken_cmd, logger, 'debug')
#                 if cluster == "cluster":
#                     cmd = cmd + "\n" + kraken_cmd
#                     krona_cmd = krona_visualization(file_prefix, Config, logger, kraken_directory, cluster)
#                     cmd = cmd + "\n" + krona_cmd
#                     generate_cluster_jobs(cmd, file_prefix, Config, logger)
#                 elif cluster == "local":
#                     call(kraken_cmd, logger)
#                     krona_cmd = krona_visualization(file_prefix, Config, logger, kraken_directory, cluster)
#                 elif cluster == "parallel-local":
#                     parallel_local_cmds.append(kraken_cmd)
#                     krona_cmd = krona_visualization(file_prefix, Config, logger, kraken_directory, cluster)
#                     parallel_local_cmds_krona.append(krona_cmd)
#     if cluster == "parallel-local":
#         if len(parallel_local_cmds) > 1:
#             complete = run_parallel(parallel_local_cmds)
#         else:
#             call(kraken_cmd, logger)
#         if len(parallel_local_cmds_krona) > 1:
#             complete = run_parallel(parallel_local_cmds_krona)
#         else:
#             call(krona_cmd, logger)
#     elif cluster == "cluster":
#         print "No support for Kraken cluster jobs. Should be run in parallel-local or local mode.\n"
#         exit()
#     elif cluster == "local":
#         #call(kraken_cmd, logger)
#         #keep_logging(kraken_cmd, kraken_cmd, logger, 'info')
#         #keep_logging(kraken_cmd, kraken_cmd, logger, 'info')
#         #krona_cmd = krona_visualization(file_prefix, Config, logger, kraken_directory, cluster)
#         #keep_logging(krona_cmd, krona_cmd, logger, 'info')
#         for i in parallel_local_cmds:
#             print i
