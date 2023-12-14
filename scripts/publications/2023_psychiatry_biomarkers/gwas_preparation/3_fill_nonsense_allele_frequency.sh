#!/bin/bash

################################################################################
# Author: T. Cameron Waller
# Date, first execution: 14 November 2023
# Date, last execution: 14 December 2023
# Date, review: 27 November 2023
################################################################################
# Note

# The purpose of this script is to fill with nonsense "0.5" allele frequency
# for any SNPs that have missing values.
# SNPs with missing values of allele frequency will otherwise be lost in the
# GWAS2VCF procedure.

# This preparation of GWAS summary statistics is for analysis in LDSC, which
# does not use information about allele frequency.

# After running this script, check the product directory to make sure that this
# procedure wrote the appropriate files again, later than all others.

################################################################################
# Organize paths.

# Identifiers or designators of parameter version and preparation batch.
identifier_preparation="tcw_2023-12-14_test"
identifier_parameter="tcw_2023-12-14_test"

# Directories.
cd ~/paths
path_directory_process=$(<"./process_psychiatric_metabolism.txt")
path_directory_dock="${path_directory_process}/dock"
path_directory_parameters="${path_directory_dock}/parameters/psychiatric_metabolism"
path_directory_source="${path_directory_dock}/gwas_preparation_${identifier_preparation}/1_gwas_format_standard"
#path_directory_source="${path_directory_dock}/gwas_preparation_${identifier_preparation}/2_gwas_assembly_grch37"
path_directory_product="${path_directory_dock}/gwas_preparation_${identifier_preparation}/3_gwas_fill_nonsense_allele_frequency"

# Files.
path_file_table_parameter="${path_directory_parameters}/table_gwas_translation_${identifier_parameter}.tsv"

# Scripts.
path_directory_partner_scripts="${path_directory_process}/partner/scripts"
path_script_process="${path_directory_partner_scripts}/gwas_clean/fill_nonsense_allele_frequency.sh"
path_script_driver="${path_directory_partner_scripts}/gwas_clean/drive_process_over_gwas.sh"

# Initialize directories.
rm -r $path_directory_product # caution
mkdir -p $path_directory_product
cd $path_directory_product

################################################################################
# Organize parameters.

report="true"

################################################################################
# Execute procedure.

/usr/bin/bash $path_script_driver \
$path_file_table_parameter \
$path_directory_source \
$path_directory_product \
$path_script_process \
$report



#
