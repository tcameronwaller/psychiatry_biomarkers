#!/bin/bash

################################################################################
# Author: T. Cameron Waller
# Date, first execution: 18 April 2022
# Date, last execution: 18 April 2023
# Review: TCW; 18 April 2023
################################################################################
# Note



################################################################################
# Organize arguments.

path_file_gwas_source=${1} # full directory path and file name for source GWAS summary statistics in standard format with GZip compression
path_file_gwas_product=${2} # full directory path and file name for product GWAS summary statistics in GCTB and GCTA-COJO ".ma" format without compression
path_directory_ld_matrix=${3} # full path to parent directory for chromosome-specific Linkage Disequilibrium (LD) reference matrices in GCTB format
name_file_ld_matrix_prefix=${4} # prefix of name of file for chromosome-specific LD reference matrices
name_file_ld_matrix_suffix=${5} # suffix of name of file for chromosome-specific LD reference matrices
path_directory_product=${6} # full path to parent directory for product files from GCTB SBayesR
name_file_effect_prefix=${7} # prefix of name of file for product files from GCTB SBayesR
name_file_effect_suffix=${8} # suffix of name of file for product files from GCTB SBayesR
observations_variant=${9} # logical binary indicator of whether counts of observations are reliable and specific to each variant (SNP)
report=${10} # whether to print reports



################################################################################
# Organize paths.

# Directories.
cd ~/paths
path_directory_tools=$(<"./waller_tools.txt")
path_gctb=$(<"./tools_gctb.txt")
path_directory_process=$(<"./process_psychiatric_metabolism.txt")

# Scripts.
path_script_gwas_format="${path_directory_process}/promiscuity/scripts/gctb/constrain_translate_gwas_standard_to_gctb.sh"
path_script_submit_batch="${path_directory_process}/promiscuity/scripts/gctb/1_submit_batch_gctb_sbayesr_chromosomes.sh"
path_script_batch_run_sbayesr="${path_directory_process}/promiscuity/scripts/gctb/2_run_batch_gctb_sbayesr.sh"
path_script_run_sbayesr="${path_directory_process}/promiscuity/scripts/gctb/run_gctb_sbayesr.sh"

# Initialize directories.
cd $path_directory_product



################################################################################
# Organize parameters.

chromosome_x="false"
threads=1

###########################################################################
# Execute procedure.

##########
# 1. Prepare GWAS summary statistics.

if true; then
  /usr/bin/bash $path_script_gwas_format \
  $path_file_gwas_source \
  $path_file_gwas_product \
  $report
fi

##########
# 2. Prepare and submit batch of jobs for processing on each chromosome.

if true; then
  /usr/bin/bash $path_script_submit_batch \
  $path_file_gwas_product \
  $path_directory_ld_matrix \
  $name_file_ld_matrix_prefix \
  $name_file_ld_matrix_suffix \
  $path_directory_product \
  $name_file_effect_prefix \
  $name_file_effect_suffix \
  $observations_variant \
  $chromosome_x \
  $path_script_batch_run_sbayesr \
  $path_script_run_sbayesr \
  $path_gctb \
  $threads \
  $report
fi



#
