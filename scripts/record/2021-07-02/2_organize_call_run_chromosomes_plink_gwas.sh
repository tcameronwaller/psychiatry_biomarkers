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
# Allow simultaneous processes for this count of GWAS studies.
# For large cohorts (20,000 - 500,000), limit to 10-20 total simultaneous GWAS
# on NCSA.
# Beyond about 15 simultaneous GWAS, PLINK2 begins to use more than 2 TB storage.
#$ -tc 13

# http://gridscheduler.sourceforge.net/htmlman/htmlman1/qsub.html

################################################################################
# Organize argument variables.

path_batch_instances=${1} # text list of information for each instance in batch
batch_instances_count=${2} # count of instances in batch
path_cohorts_models=${3} # full path to parent directory for tables of variables across cohorts
path_gwas=${4} # full path to parent directory for GWAS summary statistics
path_scripts_record=${5} # full path to parent directory for date-specific record scripts

###########################################################################
# Organize variables.

# Determine batch instance.
batch_index=$((SGE_TASK_ID-1))
readarray -t batch_instances < $path_batch_instances
instance=${batch_instances[$batch_index]}

# Separate fields from instance.
IFS=";" read -r -a array <<< "${instance}"
phenotype="${array[0]}"
cohort_model="${array[1]}"
table_cohort_model="${array[2]}"
covariates="${array[3]}"

# Organize variables.
cohort_model_phenotype="${cohort_model}_${phenotype}"
path_table_phenotypes_covariates="${path_cohorts_models}/${table_cohort_model}_${phenotype}.tsv"
path_study_gwas="${path_gwas}/${cohort_model_phenotype}"

# Report.
echo "phenotype: " ${phenotype}
echo "cohort_model: " ${cohort_model}
echo "table: " ${path_table_phenotypes_covariates}
echo "covariates: " ${covariates}

# General parameters.
threads=32 # 32, needs to match the "pe threaded" argument to scheduler
maf=0.01
chromosomes=22 # 22 # Count of chromosomes on which to run GWAS

# Echo each command to console.
#set -x
# Suppress echo each command to console.
set +x

# Initialize directories.
#rm -r $path_cohort_phenotype
if [ ! -d $path_study_gwas ]; then
    # Directory does not already exist.
    # Create directory.
    mkdir -p $path_study_gwas
fi

###########################################################################
# Execute procedure.

path_report=$path_study_gwas
analysis="${cohort_model_phenotype}"
/usr/bin/bash "${path_scripts_record}/3_run_chromosomes_plink_gwas.sh" \
$path_table_phenotypes_covariates \
$path_report \
$analysis \
$phenotype \
$covariates \
$threads \
$maf \
$chromosomes
