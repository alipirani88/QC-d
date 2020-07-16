import os
import subprocess
import re
from modules.log_modules import keep_logging
from logging_subprocess import *
from config_settings import ConfigSectionMap

def summary(filenames_array, Config, logger, prefix, output_folder):
    f = open("%s/summary.tsv" % (output_folder, prefix), 'w+')
    f.write("Sample\tRaw Coverage\tMean Coverage/Read Depth\t% bases supported by atleast 1 read\t% bases supported by atleast 5 reads\t% bases supported by atleast 9 reads\t% bases supported by atleast 10 reads\tST\tSpecies" + '\n')
    for file in filenames_array:
        filename_base = os.path.basename(file)
        if "R1_001_final.fastq.gz" in filename_base or "R1.fastq.gz" in filename_base or "1_combine.fastq.gz" in filename_base or "1_sequence.fastq.gz" in filename_base or "_forward.fastq.gz" in filename_base or "R1_001.fastq.gz" in filename_base or "_1.fastq.gz" in filename_base or ".1.fastq.gz" in filename_base or "_R1.fastq.gz" in filename_base or "_L001_R1_001.fastq.gz" in filename_base:
            # Forward reads file name and get analysis name from its name
            first_file = file
            # Get the name of reverse reads files
            if "R1_001_final.fastq.gz" in filename_base:
                second_part = filename_base.replace("R1_001_final.fastq.gz", "R2_001_final.fastq.gz")
                first_part_split = filename_base.split('R1_001_final.fastq.gz')
                first_part = first_part_split[0].replace('_L001', '')
                first_part = re.sub("_S[0-9].*_", "", first_part)
            elif "R1_001.fastq.gz" in filename_base:
                second_part = filename_base.replace("R1_001.fastq.gz", "R2_001.fastq.gz")
                first_part_split = filename_base.split('R1_001.fastq.gz')
                first_part = first_part_split[0].replace('_L001', '')
                first_part = re.sub("_S[0-9].*_", "", first_part)
            elif "_R1.fastq.gz" in filename_base:
                second_part = filename_base.replace("_R1.fastq.gz", "_R2.fastq.gz")
                first_part_split = filename_base.split('_R1.fastq.gz')
                first_part = first_part_split[0].replace('_L001', '')
                first_part = re.sub("_S[0-9].*_", "", first_part)
            elif "R1.fastq.gz" in filename_base:
                second_part = filename_base.replace("R1.fastq.gz", "R2.fastq.gz")
                first_part_split = filename_base.split('R1.fastq.gz')
                first_part = first_part_split[0].replace('_L001', '')
                first_part = re.sub("_S[0-9].*_", "", first_part)
            elif "1_combine.fastq.gz" in filename_base:
                second_part = filename_base.replace("1_combine.fastq.gz", "2_combine.fastq.gz")
                first_part_split = filename_base.split('1_combine.fastq.gz')
                first_part = first_part_split[0].replace('_L001', '')
                first_part = re.sub("_S[0-9].*_", "", first_part)
            elif "1_sequence.fastq.gz" in filename_base:
                second_part = filename_base.replace("1_sequence.fastq.gz", "2_sequence.fastq.gz")
                first_part_split = filename_base.split('1_sequence.fastq.gz')
                first_part = first_part_split[0].replace('_L001', '')
                first_part = re.sub("_S[0-9].*_", "", first_part)
            elif "_forward.fastq.gz" in filename_base:
                second_part = filename_base.replace("_forward.fastq.gz", "_reverse.fastq.gz")
                first_part_split = filename_base.split('_forward.fastq.gz')
                first_part = first_part_split[0].replace('_L001', '')
                first_part = re.sub("_S[0-9].*_", "", first_part)
            elif "R1_001.fastq.gz" in filename_base:
                second_part = filename_base.replace("R1_001.fastq.gz", "R2_001.fastq.gz")
                first_part_split = filename_base.split('R1_001.fastq.gz')
                first_part = first_part_split[0].replace('_L001', '')
                first_part = re.sub("_S[0-9].*_", "", first_part)
            elif "_1.fastq.gz" in filename_base:
                second_part = filename_base.replace("_1.fastq.gz", "_2.fastq.gz")
                first_part_split = filename_base.split('_1.fastq.gz')
                first_part = first_part_split[0].replace('_L001', '')
                first_part = re.sub("_S[0-9].*_", "", first_part)
            elif ".1.fastq.gz" in filename_base:
                second_part = filename_base.replace(".1.fastq.gz", ".2.fastq.gz")
                first_part_split = filename_base.split('.1.fastq.gz')
                first_part = first_part_split[0].replace('_L001', '')
                first_part = re.sub("_S[0-9].*_", "", first_part)


            raw_cov = sp.Popen(["grep '%s' %s/Test_Final_Coverage.txt | cut -f16" % (filename_base, output_folder)], stdout=sp.PIPE, shell=True)
            (raw_cov, err) = raw_cov.communicate()
            raw_cov = raw_cov.strip()

            cov_depth = sp.Popen(["grep '%s' %s/%s_Coverage_depth/%s_depth_of_coverage.sample_summary | cut -f3,7,8,9,10" % (filename_base, output_folder, prefix, first_part)], stdout=sp.PIPE, shell=True)
            (cov_depth, err) = cov_depth.communicate()
            cov_depth = cov_depth.strip()

            kraken_sp = sp.Popen(["awk '$4 == \"S\" {print $6,$7}' %s/Test_Kraken_results/%s*_kraken_report.txt | head -n1 | cut -f6" % (output_folder, first_part)], stdout=sp.PIPE, shell=True)
            (kraken_sp, err) = kraken_sp.communicate()
            kraken_sp = kraken_sp.strip()

            st = sp.Popen(["tail -n1 %s/%s_MLST_results/%s/mlst_report.tsv | cut -f1" % (output_folder, prefix, first_part)], stdout=sp.PIPE, shell=True)
            (st, err) = st.communicate()
            st = st.strip()

            summary_string = filename_base + "\t" + raw_cov + "\t" + cov_depth + "\t" + st + "\t" + kraken_sp + '\n'

            f.write(summary_string)