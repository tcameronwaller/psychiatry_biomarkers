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
#$ -l h_vmem=5G
# Concurrent threads; assigns value to variable NSLOTS.
# Important to specify 32 threads to avoid inconsistency with interactive
# calculations.
#$ -pe threaded 32
# Range of indices.
# Specify as argument when calling qsub.
# Array batch indices cannot start at zero.
### -t 1-100:1
##########
# Limit on concurrent processes.
# Allow simultaneous processes for this count of GWAS studies.
# For large cohorts (20,000 - 500,000), limit to 10-20 total simultaneous GWAS
# on NCSA.
# At this scale, about 15 simultaneous GWAS in PLINK2 require approximately 1.5
# Terabytes of storage (not memory).
# PLINK2 seems to require the most storage during the initial reads of variants
# such that simultaneous starts to multiple GWAS cause the requirement to peak
# temporarily.
# Also remember that as large batches complete, the final GWAS summary
# statistics accumulate and consume some of the available storage space.
#$ -tc 18

# http://gridscheduler.sourceforge.net/htmlman/htmlman1/qsub.html

################################################################################
# Organize argument variables.

path_batch_instances=${1} # text list of information for each instance in batch
batch_instances_count=${2} # count of instances in batch
path_stratification_tables=${3} # full path to parent directory for tables of phenotype and covariate variables
path_gwas_container=${4} # full path to parent directory for GWAS summary statistics
path_scripts_record=${5} # full path to parent directory for date-specific record scripts
path_process=${6} # full path to directory for all processes relevant to current project

###########################################################################
# Organize variables.

# Determine batch instance.
batch_index=$((SGE_TASK_ID-1))
readarray -t batch_instances < $path_batch_instances
instance=${batch_instances[$batch_index]}

# Separate fields from instance.
IFS=";" read -r -a array <<< "${instance}"
name_study="${array[0]}"
table_phenotypes_covariates="${array[1]}"
phenotypes="${array[2]}"
covariates="${array[3]}"

# Organize variables.
path_gwas_study_parent="${path_gwas_container}/${name_study}"
path_table_phenotypes_covariates="${path_stratification_tables}/${table_phenotypes_covariates}"

# Report.
echo "name of study: " ${name_study}
echo "path to gwas reports: " ${path_gwas_study_parent}
echo "table: " ${path_table_phenotypes_covariates}
echo "phenotypes: " ${phenotypes}
echo "covariates: " ${covariates}

# General parameters.
threads=32 # 32, needs to match the "pe threaded" argument to scheduler
maf=0.0 # run on all SNPs and filter in subsequent analyses

# Echo each command to console.
#set -x
# Suppress echo each command to console.
set +x

# Initialize directories.
rm -r $path_gwas_study_parent # be careful!
if [ ! -d $path_gwas_study_parent ]; then
    # Directory does not already exist.
    # Create directory.
    mkdir -p $path_gwas_study_parent
fi

# Scripts.
path_promiscuity_scripts="${path_process}/promiscuity/scripts"
path_scripts_gwas_process="${path_promiscuity_scripts}/gwas_process"
path_script_run_gwas="${path_scripts_gwas_process}/run_gwas_chromosomes_plink_association.sh"

###########################################################################
# Execute procedure.

/usr/bin/bash "${path_script_run_gwas}" \
$path_table_phenotypes_covariates \
$path_gwas_study_parent \
$name_study \
$phenotypes \
$covariates \
$threads \
$maf



#####
