#!/bin/bash

################################################################################
# Author: T. Cameron Waller
# Date, first execution: 30 December 2023
# Date, last execution: 30 December 2023
# Date, review: 30 December 2023
################################################################################
# Note

# count of files: 40

################################################################################
# Organize paths.

# Identifiers or designators of parameter version and preparation batch.
identifier_preparation="gwas_2023-12-22_dbsnp_rsid"
identifier_parameter="tcw_2023-12-22_dbsnp_rsid"

# Directories.
cd ~/paths
path_directory_gwas_summaries=$(<"./gwas_summaries_waller_metabolism.txt")
path_directory_process=$(<"./process_psychiatric_metabolism.txt")
path_directory_dock="${path_directory_process}/dock"

path_directory_parent_source=$path_directory_dock
path_directory_parent_product="${path_directory_gwas_summaries}/organization"

path_directory_source="${path_directory_dock}/${identifier_preparation}"
path_directory_product="${path_directory_gwas_summaries}/organization/${identifier_preparation}"

# Initialize directories.
#rm -r $path_directory_product # caution

################################################################################
# Organize parameters.

report="true"
set -x

################################################################################
# Execute procedure.

cd $path_directory_parent_product
pwd
du -hs ./*

cd $path_directory_source
pwd
du -hs ./*

cp -r $path_directory_source $path_directory_product

cd $path_directory_product
pwd
du -hs ./*

echo "----------"
echo "Count of files in batch:"
ls $path_directory_product | wc -l
echo "----------"

#rm -r $path_directory_source # caution

#
