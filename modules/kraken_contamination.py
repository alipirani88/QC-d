__author__ = 'alipirani'
import os
import subprocess
from modules.log_modules import keep_logging
from logging_subprocess import *
from config_settings import ConfigSectionMap
from modules.krona_visualization import *
from modules.generate_cluster_jobs import *
from modules.run_parallel import *


def downsample_reads(file, file2, genome_size, logger):

    # Extract basic fastq reads stats with seqtk

    gsize = int(genome_size)

    keep_logging('',
                 'Using Genome Size: %s to calculate coverage' % gsize, logger, 'debug')

    seqtk_check = "seqtk fqchk -q3 %s > /tmp/%s_fastqchk.txt" % (file, os.path.basename(file))

    keep_logging('',
                 'Running seqtk to extract Fastq statistics: %s' % seqtk_check, logger, 'debug')

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


    keep_logging('',
                 'Average Read Length: %s' % avg_len, logger, 'debug')

    keep_logging('',
                 'Total number of bases in fastq: %s' % total_bases, logger, 'debug')

    # Calculate original depth and check if it needs to be downsampled to a default coverage.
    ori_coverage_depth = int(total_bases / gsize)

    keep_logging('',
                 'Original Covarage Depth: %s x' % ori_coverage_depth, logger, 'debug')

    proc = sp.Popen(["nproc"], stdout=sp.PIPE, shell=True)
    (nproc, err) = proc.communicate()
    nproc = nproc.strip()

    if ori_coverage_depth >= 100:
        # Downsample to 100
        factor = float(100 / float(ori_coverage_depth))
        r1_sub = "/tmp/%s" % os.path.basename(file)

        # Downsample using seqtk
        try:
            keep_logging("", "Generating seqtk Downsampling command", logger, 'info')
            keep_logging("", "seqtk sample %s %s | pigz --fast -c -p %s > /tmp/%s" % (
                file, factor, nproc, os.path.basename(file)), logger, 'debug')

            seqtk_downsample = "seqtk sample %s %s | pigz --fast -c -p %s > /tmp/%s" % (
                file, factor, nproc, os.path.basename(file))
            #call(seqtk_downsample, logger)
        except sp.CalledProcessError:
            keep_logging('Error running seqtk for downsampling raw fastq reads.',
                         'Error running seqtk for downsampling raw fastq reads.', logger, 'info')
            sys.exit(1)

        if file2:
            r2_sub = "/tmp/%s" % os.path.basename(file2)
            try:
                keep_logging("", "seqtk sample %s %s | pigz --fast -c -p %s > /tmp/%s" % (
                    file2, factor, nproc, os.path.basename(file2)), logger, 'debug')
                #call("seqtk sample %s %s | pigz --fast -c -p %s > /tmp/%s" % (file2, factor, nproc, os.path.basename(file2)), logger)
            except sp.CalledProcessError:
                keep_logging('Error running seqtk for downsampling raw fastq reads.',
                             'Error running seqtk for downsampling raw fastq reads.', logger, 'exception')
                sys.exit(1)
        else:
            r2_sub = "None"

    elif ori_coverage_depth < 100:
        r1_sub = file
        r2_sub = file2
        seqtk_downsample = "cp %s /tmp/" % r1_sub


    return r1_sub, r2_sub, seqtk_downsample


def kraken_contamination(filenames_array, Config, logger, output_folder, type, samples, kraken_directory, cluster, downsample, scheduler, genome_size, dryrun):
    parallel_local_cmds = []
    parallel_local_cmds_krona = []
    cmd = ""

    for file in filenames_array:
        #file_prefix = kraken_directory + "/" + os.path.basename(file)[0:20]
	    file_prefix = kraken_directory + "/" + os.path.basename(file.replace('_R1.*fastq.gz',''))
        file2 = file.replace('_R1_', '_R2_')

        read1, read2, seqtk_downsample = downsample_reads(file, file2, genome_size, logger)
        file = read1
        keep_logging("Using downsampled reads for Kraken - %s" % file, "Using downsampled reads for Kraken - %s" % file, logger, 'info')

        kraken_cmd = seqtk_downsample
        if file.endswith('.gz'):
            file = "/tmp/%s" % os.path.basename(file)
            kraken_cmd = kraken_cmd + '\n' + "kraken --quick --fastq-input --gzip-compressed --db %s --output %s_kraken %s" % (ConfigSectionMap("kraken", Config)['db_path'],
            file_prefix, file)
            keep_logging("", kraken_cmd, logger, 'info')
            kraken_out = file_prefix + "_kraken"
            report_cmd = "kraken-report --db %s %s > %s_report.txt" % (
            ConfigSectionMap("kraken", Config)['db_path'], kraken_out, kraken_out)
            kraken_commands =  kraken_cmd + "\n" + report_cmd
            if cluster == "cluster":
                cmd = cmd + "\n" + kraken_commands
                krona_cmd = krona_visualization(file_prefix, Config, logger, kraken_directory, cluster)
                report_cmd = "kraken-report --db %s %s > %s_report.txt" % (
                ConfigSectionMap("kraken", Config)['db_path'], kraken_out, kraken_out)
                kraken_commands = kraken_commands + "\n" + krona_cmd + "\n" + report_cmd
                job_filename = generate_cluster_jobs(kraken_commands, file_prefix, scheduler, Config, logger)
                if dryrun:
                    print "Submitting job - %s\n" % job_filename
                else:
                    os.system("sbatch %s" % job_filename)
            elif cluster == "local":
                call(kraken_cmd, logger)
                call(report_cmd, logger)
                krona_cmd = krona_visualization(file_prefix, Config, logger, kraken_directory, cluster)
                call(krona_cmd, logger)
