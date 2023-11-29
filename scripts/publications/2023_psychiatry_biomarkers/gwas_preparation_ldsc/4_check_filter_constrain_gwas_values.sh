#!/bin/bash

################################################################################
# Author: T. Cameron Waller
# Date, first execution: 15 November 2023
# Date, last execution: 29 November 2023
# Date, review: 29 November 2023
################################################################################
# Note

# To redirect both standard output and standard error to the same file, use
# "&> ./standard_output_error.txt"

# $ nohup bash {file_script.sh} &> {path/file_standard_out_error.txt} &
# To monitor progress in the output, use "tail -f standard_output_error.txt".
# To stop a nohup process, use "kill -9 {process identifier}".

# Perform this procedure judiciously.
# The purpose of this procedure is to run a few final checks to prepare GWAS
# summary statistics for analysis in LDSC, LDpred2, or SBayesR.
# 1. Remove any records (rows) with empty values (cells).
# 2. Constrain probabilities (p-values) to double-float precision
# (1E-308 to 1.0).
# 3. Remove records for any SNPs with missing counts of observations, as these
# SNPs raise an error in LDpred2.

# Note: TCW; 27 November 2023
# In the filters, the following studies lost some proportion of their original
# SNPs.
# A few of these studies lost 0.1% or more of their original SNPs.
# 37872160_williams_2023
# 34226706_barton_2021_albumin
# 34002096_mullins_2021_bd_no_ukb
# 34002096_mullins_2021_bd_1
# 33587031_sinnott-armstrong_2021_testosterone_secondary_female
# 33587031_sinnott-armstrong_2021_testosterone_secondary_male
# 33587031_sinnott-armstrong_2021_testosterone_bioavailable_female
# 33587031_sinnott-armstrong_2021_testosterone_bioavailable_male
# 33587031_sinnott-armstrong_2021_shbg_female
# 33587031_sinnott-armstrong_2021_shbg_male
# 32769997_zhou_2020_thyroid_hormone
# 32099098_polimanti_2020_eur_opioid_dep_exposed
# 32099098_polimanti_2020_eur_opioid_dep_unexposed
# 32099098_polimanti_2020_eur_opioid_exposure
# 31748690_purves_2020_ukb
# 29700475_wray_2018_pgc_ukb
# 00000000_neale_2020_hypothyroidism_self
# 00000000_neale_2020_hypothyroidism_icd
# 00000000_neale_2020_hyperthyroidism_self
# 00000000_neale_2020_hyperthyroidism_icd
# 00000000_neale_2020_vitamin_d
# 00000000_neale_2020_shbg
# 00000000_neale_2020_albumin

################################################################################
# Organize paths.

# Directories.
cd ~/paths
path_directory_process=$(<"./process_psychiatric_metabolism.txt")
path_directory_dock="${path_directory_process}/dock"
path_directory_parameters="${path_directory_dock}/parameters/psychiatric_metabolism"
path_directory_source="${path_directory_dock}/gwas_preparation_ldsc_tcw_2023-11-26/3_gwas_fill_nonsense_allele_frequency"
path_directory_product="${path_directory_dock}/gwas_preparation_ldsc_tcw_2023-11-26/4_filter_constrain_gwas_values"
path_directory_log="${path_directory_product}/log"
path_directory_parent_temporary="${path_directory_process}/temporary_check_7"

# Files.
path_file_table_parameter="${path_directory_parameters}/table_gwas_translation_tcw_2023-11-26.tsv"

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
set -o xtrace

################################################################################
# Execute procedure.

if false; then
  nohup srun --chdir $path_directory_batch \
  --partition=cpu-med --nodes=1 --ntasks-per-node=1 --time=1-00:00:00 \
  --pty bash -i $path_script_driver \
  $path_file_table_parameter \
  $path_directory_source \
  $path_directory_product \
  $path_script_process \
  $report \
  1> "${path_directory_log}/standard_output.txt" \
  2> "${path_directory_log}/standard_error.txt"
fi

/usr/bin/bash $path_script_driver \
$path_file_table_parameter \
$path_directory_source \
$path_directory_product \
$path_script_process \
$report



# Remove temporary, intermediate files.
rm -r $path_directory_parent_temporary



################################################################################
# Report.
if [[ "$report" == "true" ]]; then
  echo "----------"
  echo "Script complete:"
  echo "4_check_filter_constrain_gwas_values.sh"
  echo "----------"
fi



#
