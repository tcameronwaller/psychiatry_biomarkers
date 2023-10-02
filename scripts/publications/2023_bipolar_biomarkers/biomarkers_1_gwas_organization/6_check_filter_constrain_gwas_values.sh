#!/bin/bash

################################################################################
# Author: T. Cameron Waller
# Date, first execution: 2 October 2023
# Date, last execution: 2 October 2023
# Date, review: __ October 2023
################################################################################
# Note


################################################################################
# Organize paths.

# Directories.
cd ~/paths
path_directory_process=$(<"./process_psychiatric_metabolism.txt")
path_directory_parent_temporary="${path_directory_process}/temporary_check_6"
path_directory_dock="${path_directory_process}/dock"
path_directory_parameters="${path_directory_dock}/parameters/psychiatric_metabolism"
path_directory_source="${path_directory_dock}/gwas_biomarkers_tcw_2023-09-29/5_gwas_clean_gwas2vcf"
path_directory_product="${path_directory_dock}/gwas_biomarkers_tcw_2023-09-29/6_filter_constrain_gwas_values"
# Files.
# Scripts.
path_directory_partner_scripts="${path_directory_process}/partner/scripts"
path_script_check="${path_directory_partner_scripts}/gwas_clean/check_gwas_summary_values.sh"
path_script_filter="${path_directory_partner_scripts}/gwas_clean/filter_constrain_gwas_summary_values.sh"

# Initialize directories.
rm -r $path_directory_parent_temporary # caution
mkdir -p $path_directory_parent_temporary
rm -r $path_directory_product # caution
mkdir -p $path_directory_product
cd $path_directory_product

################################################################################
# Organize parameters.

report="true"

################################################################################
# Execute procedure.

##########
# Copy the GWAS summary statistics from the previous process.
# Most sets of GWAS summary statistics do not need extra processing.
# Subsequent processes on a few studies will replace the appropriate files.
cp $path_directory_source/*.txt.gz $path_directory_product

##########
# 36376304_koskeridis_2022_c_reactive_protein

# records that do not pass checks: 3 of 6,206,409 (header line); most due to values out of range in probability (p-value)
/usr/bin/bash $path_script_check \
"${path_directory_source}/36376304_koskeridis_2022_c_reactive_protein.txt.gz" \
$path_directory_parent_temporary \
$report

##########
# 34662886_backman_2021_albumin

# records that do not pass checks: 1 of 499,614 (header line)
/usr/bin/bash $path_script_check \
"${path_directory_source}/34662886_backman_2021_albumin.txt.gz" \
$path_directory_parent_temporary \
$report

##########
# 34226706_barton_2021_albumin

# records that do not pass checks: 246,495 of 5,515,076 (header line); most due to missingness in effect and standard error
/usr/bin/bash $path_script_check \
"${path_directory_source}/34226706_barton_2021_albumin.txt.gz" \
$path_directory_parent_temporary \
$report

##########
# 32769997_zhou_2020

# records that do not pass checks: 2,993 of 22,397,081 (header line); all due to missingness in count of observations
/usr/bin/bash $path_script_check \
"${path_directory_source}/32769997_zhou_2020_thyroid_hormone.txt.gz" \
$path_directory_parent_temporary \
$report

##########
# 32581359_saevarsdottir_2020

# records that raise checks: 1 of 43,529,207 (header line)
/usr/bin/bash $path_script_check \
"${path_directory_source}/32581359_saevarsdottir_2020_thyroid_autoimmunity.txt.gz" \
$path_directory_parent_temporary \
$report

# records that raise checks: 1 of 17,412,123 (header line)
# records with allele frequency of zero: 7,106,273 (header line)
/usr/bin/bash $path_script_check \
"${path_directory_source}/32581359_saevarsdottir_2020_thyroid_autoimmunity_af_impute.txt.gz" \
$path_directory_parent_temporary \
$report




##########
# Remove files of GWAS summary statistics that are known problems.

# 32581359_saevarsdottir_2020_thyroid_autoimmunity.txt.gz
# 32242144_revez_2020_vitamin_d.txt.gz



################################################################################
# Report.
if [[ "$report" == "true" ]]; then
  echo "----------"
  echo "Script complete:"
  echo "6_check_filter_constrain_gwas_values.sh"
  echo "----------"
fi



#
