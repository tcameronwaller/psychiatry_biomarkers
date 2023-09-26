#!/bin/bash

################################################################################
# Author: T. Cameron Waller
# Date, first execution: 21 September 2023
# Date, last execution: 26 September 2023
# Date, review: 26 September 2023
################################################################################
# Note

# Perform this procedure judiciously.
# The purpose of this procedure is to run a few final checks to prepare GWAS
# summary statistics for analysis in LDSC, LDpred2, or SBayesR.
# 1. Remove any records (rows) with empty values (cells).
# 2. Constrain probabilities (p-values) to double-float precision
# (1E-308 to __).
# 3. Remove records for any SNPs with missing counts of observations, as these
# SNPs raise an error in LDpred2.


# Zhou will need a filter to remove records (rows) with missing count of observations (sample size)...
# either before or after GWAS2VCF.


################################################################################
# Organize paths.

# Directories.
cd ~/paths
path_directory_process=$(<"./process_psychiatric_metabolism.txt")
path_directory_parent_temporary="${path_directory_process}/temporary_check_4"
path_directory_dock="${path_directory_process}/dock"
path_directory_parameters="${path_directory_dock}/parameters/psychiatric_metabolism"
path_directory_source="${path_directory_dock}/gwas_biomarkers_tcw_2023-09-25/3_gwas_allele_frequency"
path_directory_product="${path_directory_dock}/gwas_biomarkers_tcw_2023-09-25/4_filter_constrain_gwas_values"
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
# 32581359_saevarsdottir_2020

# records that raise checks: 1 of 43,529,207 (header line)
/usr/bin/bash $path_script_check \
"${path_directory_source}/32581359_saevarsdottir_2020_thyroid_autoimmunity.txt.gz" \
$path_directory_parent_temporary \
$report

if true; then
  # source lines: 43,529,207
  # product lines: 43,529,207
  /usr/bin/bash $path_script_filter \
  "${path_directory_source}/32581359_saevarsdottir_2020_thyroid_autoimmunity.txt.gz" \
  "${path_directory_product}/32581359_saevarsdottir_2020_thyroid_autoimmunity.txt.gz" \
  $report
fi

# records that raise checks: 7,106,273 of 17,412,123 (header line); most due to allele frequencies of zero
/usr/bin/bash $path_script_check \
"${path_directory_source}/32581359_saevarsdottir_2020_thyroid_autoimmunity_af_impute.txt.gz" \
$path_directory_parent_temporary \
$report

if true; then
  # source lines: 17,412,123
  # product lines: 17,412,123
  /usr/bin/bash $path_script_filter \
  "${path_directory_source}/32581359_saevarsdottir_2020_thyroid_autoimmunity_af_impute.txt.gz" \
  "${path_directory_product}/32581359_saevarsdottir_2020_thyroid_autoimmunity_af_impute.txt.gz" \
  $report
fi

##########
# 32769997_zhou_2020

# records that do not pass checks: 2,993 of 22,397,081 (header line)
/usr/bin/bash $path_script_check \
"${path_directory_source}/32769997_zhou_2020_thyroid_hormone.txt.gz" \
$path_directory_parent_temporary \
$report

if true; then
  # source lines: 22,397,081
  # product lines: 22,394,089
  /usr/bin/bash $path_script_filter \
  "${path_directory_source}/32769997_zhou_2020_thyroid_hormone.txt.gz" \
  "${path_directory_product}/32769997_zhou_2020_thyroid_hormone.txt.gz" \
  $report
fi



################################################################################
# Report.
if [[ "$report" == "true" ]]; then
  echo "----------"
  echo "Script complete:"
  echo "4_check_filter_constrain_gwas_values.sh"
  echo "----------"
fi



#
