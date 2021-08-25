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
# Important to specify 32 threads to avoid inconsistency with interactive
# calculations.
#$ -pe threaded 8
# Range of indices.
# Specify as argument when calling qsub.
# Array batch indices cannot start at zero.
### -t 1-100:1
# Limit on concurrent processes.
#$ -tc 30

# http://gridscheduler.sourceforge.net/htmlman/htmlman1/qsub.html

################################################################################
# Note.

# Collection and concatenation of GWAS summary statistics runs with 4 thread
# slots ("-pe threaded 4") and 1 Gigabyte of memory ("-l h_vmem=1G"), even for
# GWAS on > 300,000 persons.

# LDSC Munge of GWAS summary statistics throws a memory error with 4 thread
# slots ("-pe threaded 4") and 1 Gigabyte of memory ("-l h_vmem=1G").

# LDSC Munge of GWAS summary statistics runs with 8 thread slots
# ("-pe threaded 8") and 1 Gigabyte of memory ("-l h_vmem=1G").

################################################################################
# Organize argument variables.

path_batch_instances=${1} # text list of information for each instance in batch
batch_instances_count=${2} # count of instances in batch
path_gwas_source_container=${3} # full path to parent directories of GWAS summary statistics for each study
path_gwas_target_container=${4} # full path to parent directories of GWAS summary statistics for each study
path_scripts_record=${5} # full path to directory of scripts for a specific analysis report date

###########################################################################
# Organize variables.

# Determine batch instance.
batch_index=$((SGE_TASK_ID-1))
readarray -t batch_instances < $path_batch_instances
study=${batch_instances[$batch_index]}

###########################################################################
# Execute procedure.

# Concatenate GWAS across chromosomes.
/usr/bin/bash "${path_scripts_record}/10_drive_gwas_concatenation_format_munge_heritability.sh" \
$study \
$path_gwas_source_container \
$path_gwas_target_container
