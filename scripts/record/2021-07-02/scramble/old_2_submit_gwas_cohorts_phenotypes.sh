#!/bin/bash

#chmod u+x script.sh

###########################################################################
###########################################################################
###########################################################################
# This script organizes directories and iteration instances then submits
# script "regress_metabolite_heritability.sh" to the Sun Grid Engine.
###########################################################################
###########################################################################
###########################################################################

################################################################################
# Organize arguments.
file_name=${1} # name of file for table with phenotypes and covariates
cohort_phenotype=${2} # unique name of cohort and phenotype
phenotype=${3} # name of table's column for phenotype, dependent variable
covariates=${4} # name of table's columns for covariates, independent variables
path_scripts_record=${5} # full path to pipeline scripts
path_cohorts=${6} # full path to parent directory for cohorts' tables
path_gwas=${7} # full path to parent directory for GWAS summary statistics

################################################################################
# Organize variables.

path_table_phenotypes_covariates="${path_cohorts}/${file_name}"
path_cohort_phenotype="$path_gwas/$cohort_phenotype"

# General parameters.
threads=32 # 32, needs to match the "pe threaded" argument to scheduler
maf=0.01
chromosomes=22 # 22 # Count of chromosomes on which to run GWAS

# Echo each command to console.
#set -x
# Suppress echo each command to console.
set +x

# Initialize directories.
rm -r $path_cohort_phenotype
mkdir -p $path_cohort_phenotype

# Submit array batch to Sun Grid Engine.
# Array batch indices cannot start at zero.
echo "----------------------------------------------------------------------"
echo "Submit array batch to Sun Grid Engine."
echo "----------------------------------------------------------------------"
qsub -t 1-${chromosomes}:1 \
-o "$path_cohort_phenotype/out.txt" -e "$path_cohort_phenotype/error.txt" \
$path_scripts_record/3_run_batch_gwas.sh \
$path_table_phenotypes_covariates \
$path_cohort_phenotype \
$cohort_phenotype \
$phenotype \
$covariates \
$threads \
$maf \
