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

path_directory_source_1="${path_directory_dock}/test_sbayesr_body_mass_tcw_2023-03-21/sbayesr_1_polygenic_scores"
path_directory_source_2="${path_directory_dock}/test_sbayesr_body_mass_tcw_2023-03-21/sbayesr_2_polygenic_scores"
path_directory_source_3="${path_directory_dock}/test_sbayesr_body_mass_tcw_2023-03-21/sbayesr_3_polygenic_scores"

path_directory_product_1="${path_directory_dock}/test_sbayesr_body_mass_tcw_2023-03-21/sbayesr_1_polygenic_scores"
path_directory_product_2="${path_directory_dock}/test_sbayesr_body_mass_tcw_2023-03-21/sbayesr_2_polygenic_scores"
path_directory_product_3="${path_directory_dock}/test_sbayesr_body_mass_tcw_2023-03-21/sbayesr_3_polygenic_scores"

# Files.
name_file_source_prefix="BMI_GIANTUKB_EUR_chromosome_"
name_file_source_suffix=".sscore"
name_file_source_not=".vars" # exclude the files that are lists of SNPs used in calculation of scores
path_file_product_1="${path_directory_product_1}/BMI_GIANTUKB_EUR_combination.tsv"
path_file_product_2="${path_directory_product_2}/BMI_GIANTUKB_EUR_combination.tsv"
path_file_product_3="${path_directory_product_3}/BMI_GIANTUKB_EUR_combination.tsv"

# Scripts.
path_script_combine_scores="${path_directory_process}/promiscuity/scripts/plink/combine_sum_polygenic_scores.sh"

# Initialize directories.
#rm -r $path_directory_product_1 # Source and product directories are the same.
#rm -r $path_directory_product_2
#rm -r $path_directory_product_3
mkdir -p $path_directory_product_1
mkdir -p $path_directory_product_2
mkdir -p $path_directory_product_3
cd $path_directory_product_1

###########################################################################
# Organize parameters.

# Combine polygenic scores across chromosomes.

# 1.
cd $path_directory_product_1
/usr/bin/bash $path_script_combine_scores \
$path_directory_source_1 \
$name_file_source_prefix \
$name_file_source_suffix \
$name_file_source_not \
$path_file_product_1

# 2.
cd $path_directory_product_2
/usr/bin/bash $path_script_combine_scores \
$path_directory_source_2 \
$name_file_source_prefix \
$name_file_source_suffix \
$name_file_source_not \
$path_file_product_2

# 3.
cd $path_directory_product_3
/usr/bin/bash $path_script_combine_scores \
$path_directory_source_3 \
$name_file_source_prefix \
$name_file_source_suffix \
$name_file_source_not \
$path_file_product_3


################################################################################
# Report.

if [[ "$report" == "true" ]]; then
  echo "----------"
  echo "Script:"
  echo "16_combine_standardize_polygenic_scores.sh"
  echo "----------"
fi



#
