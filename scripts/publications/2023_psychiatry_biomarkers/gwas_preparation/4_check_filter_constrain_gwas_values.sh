#!/bin/bash

################################################################################
# Author: T. Cameron Waller
# Date, first execution: 15 November 2023
# Date, last execution: 19 December 2023
# Date, review: 19 December 2023
################################################################################
# Note

##########
# Execution with NoHup.

# To redirect both standard output and standard error to the same file, use
# "&> ./standard_output_error.txt"
# To monitor progress in the output, use "tail -f standard_output_error.txt".
# To monitor current processes, use "top -u {user identifier}".
# To kill a single nohup process, use "kill -9 {process identifier}".
# To kill all processes by a single user, use "pkill -u {user identifier}".

# $ nohup bash {file_script.sh} &> {/path/standard_out_error.txt} &

# nohup process identifier: 221970 (TCW; 29 November 2023 14:25 Eastern Time)

##########
# Other Notes

# Perform this procedure judiciously.
# The purpose of this procedure is to run a few final checks to prepare GWAS
# summary statistics for analysis in LDSC, LDpred2, or SBayesR.
# 1. Remove any records (rows) with empty values (cells).
# 2. Constrain probabilities (p-values) to double-float precision
# (1E-308 to 1.0).
# 3. Remove records for any SNPs with missing counts of observations, as these
# SNPs raise an error in LDpred2.

# Note: TCW; 29 November 2023
# There might not be substantial differences between the "prior_1" sets of GWAS
# summary statistics that went through process 4 in the batch of 26 November
# and those in the batch of 29 November 2023.
# The main difference was that for the batch of 29 November 2023, I separated
# the filter and constraint blocks of code so that all filters and all
# constraints would apply independently to each record in the GWAS summary
# statistics.
# In the batch of 26 November 2023, there was potential that some filter or
# constraint steps would be ignored if a record had already passed a conditional
# that came previously in the sequence of "else if" conditional expressions.

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
# Organize parameters.

set -o xtrace

################################################################################
# Organize paths.

# Identifiers or designators of parameter version and preparation batch.
identifier_preparation="gwas_2023-12-19_alcohol_sex_test_ldsc_2023-12-19"
identifier_parameter="tcw_2023-12-19_alcohol_sex_test"

# Directories.
cd ~/paths
path_directory_process=$(<"./process_psychiatric_metabolism.txt")
path_directory_dock="${path_directory_process}/dock"
path_directory_parameters="${path_directory_dock}/parameters/psychiatric_metabolism"
path_directory_source="${path_directory_dock}/${identifier_preparation}/3_gwas_fill_nonsense_allele_frequency"
path_directory_product="${path_directory_dock}/${identifier_preparation}/4_filter_constrain_gwas_values"

# Files.
path_file_table_parameter="${path_directory_parameters}/table_gwas_translation_${identifier_parameter}.tsv"

# Scripts.
path_directory_partner_scripts="${path_directory_process}/partner/scripts"
path_script_check="${path_directory_partner_scripts}/gwas_clean/check_gwas_summary_values.sh"
path_file_script_process="${path_directory_partner_scripts}/gwas_clean/filter_constrain_gwas_summary_values.sh"
path_file_script_driver="${path_directory_partner_scripts}/gwas_clean/drive_process_over_gwas.sh"

# Initialize directories.
rm -r $path_directory_product # caution
mkdir -p $path_directory_product
cd $path_directory_product

################################################################################
# Organize parameters.

report="true"

################################################################################
# Execute procedure.

# I have had some difficulty getting the srun to work predictably.
#  nohup srun --chdir $path_directory_batch \
#  --partition=cpu-med --nodes=1 --ntasks-per-node=1 --time=1-00:00:00 \
#  --pty bash -i $path_file_script_driver \
#  $path_file_table_parameter \
#  $path_directory_source \
#  $path_directory_product \
#  $path_file_script_process \
#  $report \
#  1> "${path_directory_log}/standard_output.txt" \
#  2> "${path_directory_log}/standard_error.txt"

/usr/bin/bash $path_file_script_driver \
$path_file_table_parameter \
$path_directory_source \
$path_directory_product \
$path_file_script_process \
$path_directory_partner_scripts \
$report



################################################################################
# Report.
if [[ "$report" == "true" ]]; then
  echo "----------"
  echo "Script complete:"
  echo "4_check_filter_constrain_gwas_values.sh"
  echo "----------"
fi



#
