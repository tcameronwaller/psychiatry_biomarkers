#!/bin/bash

################################################################################
################################################################################
################################################################################
# Author: T. Cameron Waller
# Date, last execution: 24 January 2023
# Date, first execution: 23 December 2022
################################################################################
################################################################################
################################################################################
# Note

# TODO: write temporary files to a temporary directory other than the source directory

################################################################################
################################################################################
################################################################################



################################################################################
# Organize paths.

# Directories.
cd ~/paths
path_bgzip=$(<"./tools_waller_bgzip.txt")
path_directory_process=$(<"./process_psychiatric_metabolism.txt")
path_directory_gwas_summaries=$(<"./gwas_summaries_waller_metabolism.txt")
path_directory_parent_source="${path_directory_gwas_summaries}"
path_directory_dock="${path_directory_process}/dock"
path_directory_parameters="${path_directory_dock}/parameters/psychiatric_metabolism"
path_directory_product="${path_directory_dock}/hormone_genetics/gwas_format_standard"
# Files.
path_file_translation="${path_directory_parameters}/table_gwas_translation_tcw_2023-01-24.tsv"
# Scripts.
path_promiscuity_scripts="${path_directory_process}/promiscuity/scripts"
path_directory_script="${path_promiscuity_scripts}/gwas_format/translate_format_raw_to_standard"
path_script_drive_translations="${path_promiscuity_scripts}/utility/drive_translations_gwas_to_standard_format.sh"

# Initialize directories.
rm -r $path_directory_product
mkdir -p $path_directory_product
cd $path_directory_product

###########################################################################
# Organize parameters.

report="true"

###########################################################################
# Execute procedure.
###########################################################################

# TODO:
# 1. find all GWAS sumstats files within the previous directory (prelim) in standard format.
# 2. call a "pipeline" script within "promiscuity" that connects the various steps



/usr/bin/bash $path_script_drive_translations \
$path_file_translation \
$path_directory_parent_source \
$path_directory_script \
$path_directory_product \
$path_bgzip \
$report
