# Coverage, Quality and Contamination Analysis Pipeline

***

**The pipeline takes a file containing names of single/paired end fastq files and the type of analysis to be performed:**

**coverage:** calculates raw coverage for the given fastq samples.(need genome size; therefor its better to have a filename/samples with only one type of species)

**quality:** generates fastqc quality report of all the fastq files in the filename(It also generates a multiqc report of these fastqc results)

**screen_contamination:** runs fastq screen against the reference database(make sure you have checked the database path in fastq_screen config file as well as path to these config file has to be mentioned in the pipeline's config file)

**kraken_contamination:** Run Kraken(minikraken db only) to determine the most abundant species. Useful to determine contamination

**kraken_report:** Generate user-friendly Kraken report and krona plots from Kraken results

**coverage_depth:** Determine the depth of coverage(GATK) by mapping the reads against your choice of reference genome(check the path to reference genome in pipeline's config file)

- **optional arguments:**


```
  
  -h, --help            show this help message and exit

```

- **Required arguments:**


```

  -samples SAMPLES      Filenames of forward-paired end reads. ***One per
                        line***
                        
  -config CONFIG        Path to Config file
  
  -dir DIRECTORY        Directory of Fastq Files
  
  -analysis ANALYSIS_NAMES
                        COMMA-SEPERATED ANALYSIS_NAMES[coverage, quality,
                        screen_contamination, kraken_contamination,
                        kraken_report, coverage_depth]. Ex: -analysis coverage or -analysis coverage, quality
                        
  -o OUTPUT_FOLDER      Output Path ending with output directory name to save
                        the results
                                   
  -type TYPE            Type of analysis: SE or PE
  
  -cluster CLUSTER      Run Fastq_screen and Kraken on cluster/parallel-
                        local/local. Make Sure to check if the [CLUSTER]
                        section in config file is set up correctly.
                        
  -genome_size SIZE     Estimated Genome Size
  
  -prefix PREFIX        Prefix to use to save results files
  
  -reference REFERENCE  Reference genome to use to map against for calculating
                        the depth


```

**Note: Since Fastq screen and Kraken can be resource and time intensive, the pipeline can run individual jobs on cluster(supports PBS cluster only) or can also run multiple jobs jobs on local multiple cores.**
