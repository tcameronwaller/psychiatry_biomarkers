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

batch_index=$((SGE_TASK_ID-1))
path_batch_instances=$1

batch_index=$SGE_TASK_ID
path_table_phenotypes_covariates=$1
path_report=$2
analysis=$3
phenotypes=$4
covariates=$5
threads=$6
maf=$7

###########################################################################
# Organize paths.

# Read private, local file paths.
echo "read private file path variables and organize paths..."
cd ~/paths
path_plink2=$(<"./tools_user_plink2.txt")
path_ukb_genotype=$(<"./ukbiobank_genotype.txt")


path_temporary_gwas_format="$path_destination_parent/temporary_gwas_format.txt"

###########################################################################
# Execute procedure.

# Determine batch instance.
readarray -t batch_instances < $path_batch_instances
path_file=${batch_instances[$batch_index]}

# Extract file name.
file_name="$(basename -- $path_file)"
echo "file: " $file_name
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



# Set directory.
path_chromosome="$path_report/chromosome_$chromosome"
# Determine whether the temporary directory structure already exists.
if [ ! -d $path_chromosome ]; then
    # Directory does not already exist.
    # Create directory.
    mkdir -p $path_chromosome
fi
cd $path_chromosome
