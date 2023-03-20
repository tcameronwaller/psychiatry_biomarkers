#!/bin/bash

################################################################################
# Author: T. Cameron Waller
# Date, first execution: 20 March 2023
# Date, last execution: 20 March 2023
# Review: TCW; ___
################################################################################
# Note



################################################################################
# Organize paths.

# Directories.
cd ~/paths
path_directory_process=$(<"./process_psychiatric_metabolism.txt")
path_directory_dock="${path_directory_process}/dock" # parent directory for procedural reads and writes
path_directory_source="${path_directory_dock}/test_sbayesr_body_mass_tcw_2023-03-01/test_polygenic_scores"
path_directory_product="${path_directory_dock}/test_sbayesr_body_mass_tcw_2023-03-01/test_polygenic_scores"

# Files.
name_file_source_prefix="BMI_GIANTUKB_EUR_chromosome_"
#name_file_source_suffix=".sscore"
path_file_product="${path_directory_product}/BMI_GIANTUKB_EUR_combination.tsv"

# Scripts.
path_script_combine_standardize="${path_directory_process}/promiscuity/scripts/plink/combine_standardize_polygenic_scores.sh"

# Initialize directories.
#rm -r $path_directory_product
mkdir -p $path_directory_product
cd $path_directory_product

###########################################################################
# Organize parameters.

# Calculate polygenic scores for a single chromosome.
/usr/bin/bash $path_script_combine_standardize \
$path_directory_source \
$name_file_source_prefix \
$path_file_product



################################################################################
# Report.

if [[ "$report" == "true" ]]; then
  echo "----------"
  echo "Script:"
  echo "16_combine_standardize_polygenic_scores.sh"
  echo "----------"
fi



#
