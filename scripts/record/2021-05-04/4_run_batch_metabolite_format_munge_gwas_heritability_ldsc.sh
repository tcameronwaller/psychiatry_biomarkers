#!/bin/bash

###########################################################################
# Specify arguments for qsub command.
# Note that bash does not interpret qsub parameters, which are bash comments.
# Bash will not expand variables in qsub parameters.
# Shell.
#$ -S /bin/bash
# Name of job.
#$ -N waller_metabolite
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
#$ -pe threaded 32
# Range of indices.
# Specify as argument when calling qsub.
# Array batch indices cannot start at zero.
### -t 1-100:1
# Limit on concurrent processes.
#$ -tc 50

# http://gridscheduler.sourceforge.net/htmlman/htmlman1/qsub.html

################################################################################
# Organize argument variables.

path_batch_instances=${1} # text list of information for each instance in batch
batch_instances_count=${2} # count of instances in batch
metabolite_study=${3}
metabolite_file_prefix=${4}
metabolite_file_suffix=${5}
path_genetic_reference=${6}
path_gwas_study=${7}
path_heritability_study=${8}
path_script_gwas_format=${9}
path_format_munge_gwas_heritability_ldsc=${10}
path_promiscuity_scripts=${11}
path_ldsc=${12} # full path to executable LDSC
report=${13}

###########################################################################
# Organize variables.

# Determine batch instance.
batch_index=$((SGE_TASK_ID-1))
readarray -t batch_instances < $path_batch_instances
path_source_file=${batch_instances[$batch_index]}

# Determine file name.
#file_name=$path_source_file
file_name="$(basename -- $path_source_file)"
# Determine metabolite identifier.
# Refer to documnetation for test: https://www.freebsd.org/cgi/man.cgi?test
# Bash script more or less ignores empty string argument, so it does not
# work well to pass an empty string as an argument.
# Instead use a non-empty string, such as "null".
# if [[ ! -z "$name_prefix" ]]; then
metabolite=${file_name}
if [[ "$metabolite_file_prefix" != "null" ]]; then
  metabolite=${metabolite/$metabolite_file_prefix/""}
fi
if [[ "$metabolite_file_suffix" != "null" ]]; then
  metabolite=${metabolite/$metabolite_file_suffix/""}
fi
# Report.
if [[ "$report" == "true" ]]; then
  echo "----------------------------------------------------------------------"
  echo "----------------------------------------------------------------------"
  echo "----------------------------------------------------------------------"
  echo "metabolite study: " $metabolite_study
  echo "path to metabolite file: " $path_source_file
  echo "file: " $file_name
  echo "metabolite: " $metabolite
  echo "----------"
fi

study=${metabolite_study}
name_prefix=${metabolite} # file name prefix or "null"
/usr/bin/bash "$path_format_munge_gwas_heritability_ldsc" \
$study \
$name_prefix \
$path_source_file \
$path_genetic_reference \
$path_gwas_study \
$path_heritability_study \
$path_script_gwas_format \
$path_promiscuity_scripts \
$path_ldsc \
$report
