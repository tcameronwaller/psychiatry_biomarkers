#!/bin/bash

################################################################################
# Author: T. Cameron Waller
# Date, first execution: 23 December 2022
# Date, last execution: 26 September 2023
# Date, review: 19 September 2023
################################################################################
# Note



################################################################################



################################################################################
# Organize paths.

# Directories.
cd ~/paths
path_directory_process=$(<"./process_psychiatric_metabolism.txt")
path_directory_dock="${path_directory_process}/dock"
path_directory_parameters="${path_directory_dock}/parameters/psychiatric_metabolism"

path_directory_source="${path_directory_dock}/gwas_biomarkers_tcw_2023-09-25/4_filter_constrain_gwas_values"
path_directory_product="${path_directory_dock}/gwas_biomarkers_tcw_2023-09-25/5_gwas_clean_gwas2vcf"
# Files.
path_file_translation="${path_directory_parameters}/table_gwas_translation_tcw_2023-09-19_biomarkers.tsv"
# Scripts.
path_directory_partner_scripts="${path_directory_process}/partner/scripts"
path_script_submit_batch="${path_directory_partner_scripts}/gwas_clean/1_submit_batch_pipe_gwas_clean.sh"

# Initialize directories.
rm -r $path_directory_product # caution
mkdir -p $path_directory_product
cd $path_directory_product

###########################################################################
# Organize parameters.

report="true"

###########################################################################
# Execute procedure.
###########################################################################

/usr/bin/bash $path_script_submit_batch \
$path_file_translation \
$path_directory_source \
$path_directory_product \
$report



#
