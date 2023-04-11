#!/bin/bash

################################################################################
# Author: T. Cameron Waller
# Date, first execution: 21 March 2023
# Date, last execution: 11 April 2023
# Review: TCW; 11 April 2023
################################################################################
# Note

# Source Format
# Description: Format polygenic scores standard for collection and standardization.
# File suffix: ".txt"
# File type: text
# File compression: none
# Delimiter: Tab
# Columns: identifier   score (TCW; 2023-03-21)
#          1            2

################################################################################
# Organize paths.

# Directories.
cd ~/paths
path_directory_process=$(<"./process_psychiatric_metabolism.txt")
path_directory_dock="${path_directory_process}/dock" # parent directory for procedural reads and writes

path_directory_source="${path_directory_dock}/test_sbayesr_body_mass_tcw_2023-03-21/polygenic_scores_comparison/BMI_GIANTUKB_EUR_standard"
path_directory_product="${path_directory_dock}/test_sbayesr_body_mass_tcw_2023-03-21/polygenic_scores_comparison/BMI_GIANTUKB_EUR_collection"

# Files.
name_file_source_prefix="score_" # Must not be empty string.
name_file_source_suffix=".txt" # Must not be empty string.
name_file_source_not=".place_holder" # exclude any files that include this character string in file name
path_file_product="${path_directory_product}/table_scores_collection.tsv"

# Scripts.
path_script_collect_standardize="${path_directory_process}/promiscuity/scripts/plink/collect_standardize_polygenic_scores.sh"


# Initialize directories.
rm -r $path_directory_product
mkdir -p $path_directory_product
cd $path_directory_product

################################################################################
# Organize parameters.



################################################################################
# Execute procedure.

# Collect and standardize polygenic scores.

/usr/bin/bash $path_script_collect_standardize \
$path_directory_source \
$name_file_source_prefix \
$name_file_source_suffix \
$name_file_source_not \
$path_file_product




#
