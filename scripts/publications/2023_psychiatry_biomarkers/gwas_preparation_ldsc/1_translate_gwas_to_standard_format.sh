#!/bin/bash

################################################################################
# Author: T. Cameron Waller
# Date, first execution: 23 December 2022
# Date, last execution: 13 November 2023
# Date, review: 13 November 2023
################################################################################
# Note

# count of files: __
# After format translation (and maybe before) the following study was an outlier
# for its small size in terms of disk memory usage.

# 34662886_backman_2021_albumin

################################################################################
# Organize paths.

# Directories.
cd ~/paths
path_bgzip=$(<"./tools_bgzip.txt")
path_directory_process=$(<"./process_psychiatric_metabolism.txt")
path_directory_gwas_summaries=$(<"./gwas_summaries_waller_metabolism.txt")
path_directory_parent_source="${path_directory_gwas_summaries}"
path_directory_dock="${path_directory_process}/dock"
path_directory_parameters="${path_directory_dock}/parameters/psychiatric_metabolism"
path_directory_product="${path_directory_dock}/gwas_preparation_ldsc_tcw_2023-11-13/1_gwas_format_standard"
# Files.
path_file_translation="${path_directory_parameters}/table_gwas_translation_tcw_2023-11-13.tsv"

# Scripts.
path_directory_partner_scripts="${path_directory_process}/partner/scripts"
path_directory_script="${path_directory_partner_scripts}/gwas_format/translate_format_raw_to_standard"
path_script_drive_translations="${path_directory_partner_scripts}/gwas_format/drive_translations_gwas_to_standard_format.sh"

# Initialize directories.
rm -r $path_directory_product # caution
mkdir -p $path_directory_product
cd $path_directory_product

################################################################################
# Organize parameters.

report="true"

################################################################################
# Execute procedure.

/usr/bin/bash $path_script_drive_translations \
$path_file_translation \
$path_directory_parent_source \
$path_directory_script \
$path_directory_product \
$path_bgzip \
$report



#
