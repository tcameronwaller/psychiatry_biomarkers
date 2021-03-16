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
# Priority -1023 to 1024.
# Increasing priority from default (zero) requires operator privileges.
### -p 1024
# Memory per iteration.
# Segmentation errors commonly indicate a memory error.
#$ -l h_vmem=5G
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
path_genetic_reference=$4 # full path to genetic reference access
name_prefix=$5 # file name prefix before metabolite identifier or empty string
name_suffix=$6 # file name suffix after metabolite identifier or empty string
path_script_gwas_organization=$7 # full path to script to use for format organization
path_scripts=$8 # full path to scripts for current implementation pipeline
path_promiscuity_scripts=$9 # full path to scripts from promiscuity package

# Determine batch instance.
batch_index=$((SGE_TASK_ID-1))
readarray -t batch_instances < $path_batch_instances
path_file=${batch_instances[$batch_index]}

###########################################################################
# Execute procedure.

/usr/bin/bash "$path_scripts/6_execute_procedure_metabolite.sh" \
$path_file \
$path_destination_parent \
$path_genetic_reference \
$name_prefix \
$name_suffix \
$path_script_gwas_organization \
$path_scripts \
$path_promiscuity_scripts
