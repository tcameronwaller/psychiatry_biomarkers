#!/bin/bash

###########################################################################
# Specify arguments for qsub command.
# Note that bash does not interpret qsub parameters, which are bash comments.
# Bash will not expand variables in qsub parameters.
# Shell.
#$ -S /bin/bash
# Name of job.
#$ -N waller_gwas
# Contact.
# "b": beginning, "e": end, "a": abortion, "s": suspension, "n": never
#$ -M tcameronwaller@gmail.com
#$ -m as
# Standard output and error.
# Specify as arguments when calling qsub.
### -o "./out"
### -e "./error"
# Queue.
# "1-hour", "1-day", "4-day", "7-day", "30-day", "lg-mem"
#$ -q 1-hour
# Priority 0-15.
### -p -10
# Memory per iteration.
# Segmentation errors commonly indicate a memory error.
#$ -l h_vmem=1G
# Concurrent threads; assigns value to variable NSLOTS.
#$ -pe threaded 8
# Range of indices.
# Specify as argument when calling qsub.
# Array batch indices cannot start at zero.
### -t 1-100:1
# Limit on concurrent processes.
#$ -tc 50

# http://gridscheduler.sourceforge.net/htmlman/htmlman1/qsub.html

################################################################################
# Note.

# Collection and concatenation of GWAS summary statistics runs with 4 thread
# slots ("-pe threaded 4") and 1 Gigabyte of memory ("-l h_vmem=1G"), even for
# GWAS on > 300,000 persons.

################################################################################
# Organize argument variables.

path_batch_instances=${1} # text list of information for each instance in batch
batch_instances_count=${2} # count of instances in batch
pattern_gwas_report_file=${3} # string glob pattern by which to recognize PLINK2 GWAS report files
path_gwas_source_container=${4} # full path to parent directories of GWAS summary statistics for each study
path_gwas_target_container=${5} # full path to parent directories of GWAS summary statistics for each study
path_scripts_record=${6} # full path to directory of scripts for a specific analysis report date
path_process=${7} # full path to directory for all processes relevant to current project
chromosome_x=${8} # whether to collect GWAS summary statistics report for Chromosome X

###########################################################################
# Organize variables.

# Determine batch instance.
batch_index=$((SGE_TASK_ID-1))
readarray -t batch_instances < $path_batch_instances
study=${batch_instances[$batch_index]}

###########################################################################
# Execute procedure.

# Concatenate GWAS across chromosomes.
/usr/bin/bash "${path_scripts_record}/8_drive_gwas_concatenation.sh" \
$study \
$pattern_gwas_report_file \
$path_gwas_source_container \
$path_gwas_target_container \
$path_scripts_record \
$path_process \
$chromosome_x



#
