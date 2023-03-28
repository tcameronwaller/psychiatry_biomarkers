#!/bin/bash

################################################################################
################################################################################
################################################################################
# Author: T. Cameron Waller
# Date, first execution: 27 March 2023
# Date, last execution: 27 March 2023
################################################################################
################################################################################
################################################################################
# Note



################################################################################
################################################################################
################################################################################



################################################################################
# Organize paths.

# Directories.
cd ~/paths
path_directory_process=$(<"./process_psychiatric_metabolism.txt")
path_directory_dock="${path_directory_process}/dock"
path_directory_source="${path_directory_dock}/hormone_genetics_tcw_2023-02-24/gwas_effective_observations"
path_directory_product="${path_directory_dock}/hormone_genetics_tcw_2023-02-24/plots_qq_gwas"

# Files.
name_file_source_prefix="_20" # must not be empty string
name_file_source_suffix=".txt.gz" # must not be empty string
name_file_source_not=".place_holder" # exclude any files that include this character string in file name

# Scripts.
path_file_script_create_plots="${path_directory_process}/promiscuity/scripts/gwas_clean/create_gwas_qq_plots.sh"

# Initialize directories.
rm -r $path_directory_product
mkdir -p $path_directory_product

################################################################################
# Organize parameters.

# Report.
report="true"

################################################################################
# Execute procedure.

# Collect and standardize polygenic scores.

/usr/bin/bash $path_file_script_create_plots \
$path_directory_source \
$name_file_source_prefix \
$name_file_source_suffix \
$name_file_source_not \
$path_directory_product


################################################################################
# Report.
if [[ "$report" == "true" ]]; then
  echo "----------"
  echo "Script complete:"
  echo "5_create_qq_plots.sh"
  echo "----------"
fi



#
