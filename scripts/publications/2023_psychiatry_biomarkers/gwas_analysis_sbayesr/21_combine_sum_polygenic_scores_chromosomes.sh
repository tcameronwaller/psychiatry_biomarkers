#!/bin/bash

################################################################################
# Author: T. Cameron Waller
# Date, first execution: 20 March 2023
# Date, last execution: 19 April 2023
# Review: TCW; 19 April 2023
################################################################################
# Note



################################################################################
# Organize paths.

# Directories.
cd ~/paths
path_directory_process=$(<"./process_psychiatric_metabolism.txt")
path_directory_dock="${path_directory_process}/dock" # parent directory for procedural reads and writes

path_directory_source="${path_directory_dock}/hormone_genetics_tcw_2023-02-24/sbayesr_plink_polygenic_scores_16g_per_cpu_index_test"
path_directory_product="${path_directory_dock}/hormone_genetics_tcw_2023-02-24/sbayesr_plink_polygenic_scores_combination"

# Files.
name_file_source_prefix="32769997_zhou_2020_tsh_chromosome_"
name_file_source_suffix="_tcw.sscore"
name_file_source_not=".vars" # exclude the files that are lists of SNPs used in calculation of scores
path_file_product="${path_directory_product}/score_32769997_zhou_2020_tsh.txt"

# Scripts.
path_script_combine_scores="${path_directory_process}/promiscuity/scripts/plink/combine_sum_polygenic_scores.sh"

# Initialize directories.
rm -r $path_directory_product
mkdir -p $path_directory_product
cd $path_directory_product

################################################################################
# Organize parameters.



###########################################################################
# Execute procedure.

# Combine polygenic scores across chromosomes.

/usr/bin/bash $path_script_combine_scores \
$path_directory_source \
$name_file_source_prefix \
$name_file_source_suffix \
$name_file_source_not \
$path_file_product



################################################################################
# Report.

if [[ "$report" == "true" ]]; then
  echo "----------"
  echo "Script completion:"
  echo $0 # Print full file path to script.
  echo "21_combine_sum_polygenic_scores_chromosomes.sh"
  echo "----------"
fi



#
