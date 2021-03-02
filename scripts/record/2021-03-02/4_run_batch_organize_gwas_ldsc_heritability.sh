#!/bin/bash

###########################################################################
# Specify arguments for qsub command.
# Note that bash does not interpret qsub parameters, which are bash comments.
# Bash will not expand variables in qsub parameters.
# Shell.
#$ -S /bin/bash
# Name of job.
#$ -N waller_metabolites
# Contact.
# "b": beginning, "e": end, "a": abortion, "s": suspension, "n": never
#$ -M tcameronwaller@gmail.com
#$ -m ase
# Standard output and error.
# Specify as arguments when calling qsub.
### -o "./out"
### -e "./error"
# Queue.
# "1-hour", "1-day", "4-day", "7-day", "30-day", "lg-mem"
#$ -q 1-day
# Priority 0-15.
### -p -10
# Memory per iteration.
# Segmentation errors commonly indicate a memory error.
#$ -l h_vmem=10G
# Concurrent threads; assigns value to variable NSLOTS.
# Important to specify 32 threads to avoid inconsistency with interactive
# calculations.
#$ -pe threaded 32
# Range of indices.
# Specify as argument when calling qsub.
# Array batch indices cannot start at zero.
### -t 1-100:1
# Limit on concurrent processes.
#$ -tc 30

# http://gridscheduler.sourceforge.net/htmlman/htmlman1/qsub.html

###########################################################################
# Organize variables.

path_batch_instances=$1 # text list of information for each instance in batch
batch_instances_count=$2 # count of instances in batch
path_destination_parent=$3 # full path to destination directory
name_prefix=$4 # file name prefix before metabolite identifier or empty string
name_suffix=$5 # file name suffix after metabolite identifier or empty string
path_script_gwas_organization=$6 # full path to script to use for format organization
path_scripts=$7 # full path to scripts for current implementation pipeline
path_promiscuity_scripts=$8 # full path to scripts from promiscuity package

# Determine batch instance.
batch_index=$((SGE_TASK_ID-1))
readarray -t batch_instances < $path_batch_instances
path_file=${batch_instances[$batch_index]}

# Determine file name.
file_name="$(basename -- $path_file)"
echo "file: " $file_name

# Determine metabolite identifier.
# Refer to documnetation for test: https://www.freebsd.org/cgi/man.cgi?test
metabolite=${file_name}
if [[ ! -z "$name_prefix" ]]; then
  metabolite=${metabolite/$name_prefix/""}
fi
if [[ ! -z "$name_suffix" ]]; then
  metabolite=${metabolite/$name_suffix/""}
fi


###########################################################################
# Organize paths.
# Read private, local file paths.
cd ~/paths
path_ldsc=$(<"./tools_ldsc.txt")
path_temporary_gwas_format="${path_destination_parent}/temporary_gwas_format.txt"
path_temporary_gwas_munge="${path_destination_parent}/temporary_gwas_munge.txt"

###########################################################################
# Execute procedure.

# Organize information in format for LDSC.
# Parameters.
report="true" # "true" or "false"
/usr/bin/bash "$path_scripts/3_organize_gwas_ldsc_33437055_panyard_2021.sh" \
$file_name \
$path_file \
$path_temporary_gwas_format \
$path_parent \
$path_calculate_z_score_column_5_of_6 \
$report
# Munge.
# Heritability.
