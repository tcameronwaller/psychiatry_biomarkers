#!/bin/bash

################################################################################
# Author: T. Cameron Waller
# Date, first execution: 15 November 2023
# Date, last execution: 16 November 2023
# Date, review: 16 November 2023
################################################################################
# Note

# TODO:
# Revert back to the full filter script.



# Perform this procedure judiciously.
# The purpose of this procedure is to run a few final checks to prepare GWAS
# summary statistics for analysis in LDSC, LDpred2, or SBayesR.
# 1. Remove any records (rows) with empty values (cells).
# 2. Constrain probabilities (p-values) to double-float precision
# (1E-308 to 1.0).
# 3. Remove records for any SNPs with missing counts of observations, as these
# SNPs raise an error in LDpred2.


# TODO: TCW; 16 November 2023
# I think that the filter on the allele designations might be too stringent.
# GWAS2VCF might fill in the allele designations from the rsID.
# LDSC does not actually need the allele designations.


################################################################################
# Organize paths.

# Directories.
cd ~/paths
path_directory_process=$(<"./process_psychiatric_metabolism.txt")
path_directory_dock="${path_directory_process}/dock"
path_directory_parameters="${path_directory_dock}/parameters/psychiatric_metabolism"
path_directory_source="${path_directory_dock}/gwas_preparation_ldsc_tcw_2023-11-13/3_gwas_fill_nonsense_allele_frequency"
path_directory_product="${path_directory_dock}/gwas_preparation_ldsc_tcw_2023-11-13/4_filter_constrain_gwas_values"
path_directory_parent_temporary="${path_directory_process}/temporary_check_7"

# Files.
#path_file_translation="${path_directory_parameters}/table_gwas_translation_tcw_2023-11-13.tsv"
path_file_translation="${path_directory_parameters}/table_gwas_translation_tcw_2023-11-13_alcohol_sex_hormones.tsv"

# Scripts.
path_directory_partner_scripts="${path_directory_process}/partner/scripts"
path_script_check="${path_directory_partner_scripts}/gwas_clean/check_gwas_summary_values.sh"
path_script_process="${path_directory_partner_scripts}/gwas_clean/filter_constrain_gwas_summary_values.sh"
path_script_driver="${path_directory_partner_scripts}/gwas_clean/drive_process_over_gwas.sh"

# Initialize directories.
rm -r $path_directory_product # caution
mkdir -p $path_directory_product
#mkdir -p $path_directory_parent_temporary
cd $path_directory_product

################################################################################
# Organize parameters.

report="true"

################################################################################
# Execute procedure.

/usr/bin/bash $path_script_driver \
$path_file_translation \
$path_directory_source \
$path_directory_product \
$path_script_process \
$report

# Remove temporary, intermediate files.
#rm -r $path_directory_parent_temporary



################################################################################
# Report.
if [[ "$report" == "true" ]]; then
  echo "----------"
  echo "Script complete:"
  echo "4_check_filter_constrain_gwas_values.sh"
  echo "----------"
fi



#
