#!/bin/bash

###########################################################################
# Organize paths.
# Read private, local file paths.
#echo "read private file path variables and organize paths..."
cd ~/paths
path_bgzip=$(<"./tools_waller_bgzip.txt")
path_directory_process=$(<"./process_psychiatric_metabolism.txt")
path_directory_gwas_summaries=$(<"./gwas_summaries_waller_metabolism.txt")
path_directory_parent_source="${path_directory_gwas_summaries}"
path_dock="${path_directory_process}/dock"
path_directory_parameters="${path_dock}/parameters/psychiatric_metabolism"
path_directory_product="${path_dock}/hormone_genetics/gwas_format_standard"
# Scripts.
path_promiscuity_scripts="${path_directory_process}/promiscuity/scripts"
path_directory_script="${path_promiscuity_scripts}/gwas_process/format_gwas_team"
path_script_drive_translations="${path_promiscuity_scripts}/utility/drive_translations_gwas_to_standard_format.sh"
# Parameters.
path_file_translation="${path_directory_parameters}/table_gwas_translation_tcw_2022-12-23.tsv"

# Initialize directories.
rm -r $path_directory_product
mkdir -p $path_directory_product

###########################################################################
# Organize parameters.

report="true"

###########################################################################
# Execute procedure.
###########################################################################

/usr/bin/bash $path_script_drive_translations \
$path_file_translation \
$path_directory_parent_source \
$path_directory_script \
$path_directory_product \
$path_bgzip \
$report
