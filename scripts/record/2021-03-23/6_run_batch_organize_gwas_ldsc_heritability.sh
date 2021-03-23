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
#$ -l h_vmem=3G
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

################################################################################
# Organize argument variables.

path_batch_instances=${1} # text list of information for each instance in batch
batch_instances_count=${2} # count of instances in batch
phenotype_study=${3} # identifier of GWAS study for phenotype
metabolite_study=${4} # identifier of GWAS study for metabolites
name_prefix=${5} # file name prefix before metabolite identifier or "null"
name_suffix=${6} # file name suffix after metabolite identifier or "null"
path_genetic_reference=${7} # full path to parent directory with genetic reference files for LDSC
path_phenotype_gwas=${8} # full path to parent directory for formatted GWAS summary statistics for phenotype
path_study_gwas=${9} # full path to parent directory for formatted GWAS summary statistics for metabolites in study
path_study_heritability=${10} # full path to parent directory for LDSC heritability estimation for metabolites in study
path_study_genetic_correlation=${11} # full path to parent directory for LDSC genetic correlation for metabolites in study
path_scripts_record=${12} # full path to pipeline scripts
path_script_gwas_organization=${13} # full path to script to use to organize format of GWAS summary statistics for metabolites in study
path_promiscuity_scripts=${14} # complete path to directory of scripts for z-score standardization
report=${15} # whether to print reports

###########################################################################
# Organize variables.

# Determine batch instance.
batch_index=$((SGE_TASK_ID-1))
readarray -t batch_instances < $path_batch_instances
path_source_file=${batch_instances[$batch_index]}
#source_file="$(basename -- $path_source_file)"

###########################################################################
# Execute procedure.

/usr/bin/bash "${path_scripts_record}/7_execute_procedure_metabolite.sh" \
$phenotype_study \
$metabolite_study \
$path_source_file \
$name_prefix \
$name_suffix \
$path_genetic_reference \
$path_phenotype_gwas \
$path_study_gwas \
$path_study_heritability \
$path_study_genetic_correlation \
$path_script_gwas_organization \
$path_promiscuity_scripts \
$report
