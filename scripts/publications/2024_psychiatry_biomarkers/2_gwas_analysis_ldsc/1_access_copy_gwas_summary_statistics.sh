#!/bin/bash

################################################################################
# Author: T. Cameron Waller
# Date, first execution: 25 May 2023
# Date, last execution: 8 January 2024
# Date, review: 8 January 2023
################################################################################
# Note



################################################################################
# Organize paths.

# Identifiers or designators of parameter version, preparation batch, and
# analysis batch.
identifier_analysis="gwas_2023-12-30_ldsc_2024-01-08"
identifier_parameter="tcw_2023-12-30_dbsnp_rsid"

# Directories.
cd ~/paths
path_directory_gwas_summaries=$(<"./gwas_summaries_waller_metabolism.txt")
path_directory_gwas="${path_directory_gwas_summaries}/organization"
path_directory_process=$(<"./process_psychiatric_metabolism.txt")
path_directory_dock="${path_directory_process}/dock"

name_source_directory="8_gwas_effective_observations" # Name of directory with source GWAS summary statistics.
path_directory_source="${path_directory_gwas}/gwas_preparation_2023-12-30/${name_source_directory}"
path_directory_product_parent="${path_directory_dock}/${identifier_analysis}"
path_directory_product_child="${path_directory_product_parent}/${name_source_directory}"
path_directory_product="${path_directory_product_parent}/1_gwas_summaries_source"

# Initialize directories.
rm -r $path_directory_product_parent # caution
rm -r $path_directory_product # caution
mkdir -p $path_directory_product_parent
cd $path_directory_product_parent

################################################################################
# Organize parameters.

report="true"
set -x

################################################################################
# Execute procedure.

cp -r $path_directory_source $path_directory_product_parent
mv $path_directory_product_child $path_directory_product

ls $path_directory_product | wc -l

################################################################################
# Report.
if [[ "$report" == "true" ]]; then
  echo "----------"
  echo "Script complete:"
  echo $0 # Print full file path to script.
  echo "1_access_copy_gwas_summary_statistics.sh"
  echo "----------"
fi
