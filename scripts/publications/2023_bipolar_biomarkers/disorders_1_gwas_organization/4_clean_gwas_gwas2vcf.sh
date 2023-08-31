#!/bin/bash

################################################################################
# Author: T. Cameron Waller
# Date, first execution: 23 December 2022
# Date, last execution: 31 August 2023
# Date, review: 31 August 2023
################################################################################
# Note



################################################################################
# Organize paths.

# Directories.
cd ~/paths
path_directory_process=$(<"./process_psychiatric_metabolism.txt")
path_directory_dock="${path_directory_process}/dock"
path_directory_parameters="${path_directory_dock}/parameters/psychiatric_metabolism"

path_directory_source="${path_directory_dock}/gwas_disorders_tcw_2023-08-31/3_gwas_allele_frequency"
path_directory_product="${path_directory_dock}/gwas_disorders_tcw_2023-08-31/4_gwas_clean_gwas2vcf"
# Files.
path_file_translation="${path_directory_parameters}/table_gwas_translation_tcw_2023-08-31.tsv"
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
