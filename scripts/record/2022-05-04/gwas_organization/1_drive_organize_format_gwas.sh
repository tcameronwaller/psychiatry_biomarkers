#!/bin/bash

###########################################################################
###########################################################################
###########################################################################
# ...
###########################################################################
###########################################################################
###########################################################################

# TODO: TCW; 03 June 2022
# TODO: I could parse information from original directory names and use that information to define new file names automatically...
# TODO: I will eventually want to submit a batch job for all of the GWAS summary statistics files...

################################################################################
# Organize paths.
# Read private, local file paths.
cd ~/paths
path_gwas_summaries=$(<"./gwas_summaries_waller_metabolism.txt")
path_process=$(<"./process_psychiatric_metabolism.txt")
path_dock="$path_process/dock"

path_directory_gwas_concatenation="${path_dock}/gwas_concatenation"
path_directory_gwas_format="${path_dock}/gwas_format_team"

# Scripts.
path_promiscuity_scripts="${path_process}/promiscuity/scripts"
path_script_format_gwas_linear="${path_promiscuity_scripts}/gwas_process/format_gwas_team/format_gwas_team_plink_linear.sh"
path_script_format_gwas_logistic="${path_promiscuity_scripts}/gwas_process/format_gwas_team/format_gwas_team_plink_logistic.sh"

# Initialize directories.
rm -r $path_directory_gwas_format
mkdir -p $path_directory_gwas_format

###########################################################################
# Execute procedure.

# Organize multi-dimensional array of information about studies.
studies=()
# [regression type] ; [full path to source file with GWAS summary statistics] ; [full path to product file for GWAS summary statistics]
# Sex Steroid Hormone Binding Globulin (SHBG).
studies+=(
  "linear;\
  ${path_directory_gwas_concatenation}/steroid_globulin_linear_2/female_male_priority_male_joint_1_steroid_globulin_imputation_log/gwas_concatenation.txt.gz;\
  ${path_directory_gwas_format}/ukbiobank_tcw_2022-05-04/steroid_globulin_imputation_log_female_male_joint_1.txt.gz"
)
studies+=(
  "linear;\
  ${path_directory_gwas_concatenation}/steroid_globulin_linear_2/female_male_priority_male_unadjust_steroid_globulin_imputation_log/gwas_concatenation.txt.gz;\
  ${path_directory_gwas_format}/ukbiobank_tcw_2022-05-04/steroid_globulin_imputation_log_female_male_unadjust.txt.gz"
)
studies+=(
  "linear;\
  ${path_directory_gwas_concatenation}/steroid_globulin_linear_1/female_joint_1_steroid_globulin_imputation_log/gwas_concatenation.txt.gz;\
  ${path_directory_gwas_format}/ukbiobank_tcw_2022-05-04/steroid_globulin_imputation_log_female_joint_1.txt.gz"
)
studies+=(
  "linear;\
  ${path_directory_gwas_concatenation}/steroid_globulin_linear_1/female_unadjust_steroid_globulin_imputation_log/gwas_concatenation.txt.gz;\
  ${path_directory_gwas_format}/ukbiobank_tcw_2022-05-04/steroid_globulin_imputation_log_female_unadjust.txt.gz"
)
studies+=(
  "linear;\
  ${path_directory_gwas_concatenation}/steroid_globulin_linear_1/male_joint_1_steroid_globulin_imputation_log/gwas_concatenation.txt.gz;\
  ${path_directory_gwas_format}/ukbiobank_tcw_2022-05-04/steroid_globulin_imputation_log_male_joint_1.txt.gz"
)
studies+=(
  "linear;\
  ${path_directory_gwas_concatenation}/steroid_globulin_linear_1/male_unadjust_steroid_globulin_imputation_log/gwas_concatenation.txt.gz;\
  ${path_directory_gwas_format}/ukbiobank_tcw_2022-05-04/steroid_globulin_imputation_log_male_unadjust.txt.gz"
)
# Testosterone.
studies+=(
  "linear;\
  ${path_directory_gwas_concatenation}/testosterone_linear/female_joint_1_testosterone_imputation_log/gwas_concatenation.txt.gz;\
  ${path_directory_gwas_format}/ukbiobank_tcw_2022-05-04/testosterone_imputation_log_female_joint_1.txt.gz"
)
studies+=(
  "linear;\
  ${path_directory_gwas_concatenation}/testosterone_linear/female_unadjust_testosterone_imputation_log/gwas_concatenation.txt.gz;\
  ${path_directory_gwas_format}/ukbiobank_tcw_2022-05-04/testosterone_imputation_log_female_unadjust.txt.gz"
)
studies+=(
  "linear;\
  ${path_directory_gwas_concatenation}/testosterone_linear/male_joint_1_testosterone_imputation_log/gwas_concatenation.txt.gz;\
  ${path_directory_gwas_format}/ukbiobank_tcw_2022-05-04/testosterone_imputation_log_male_joint_1.txt.gz"
)
studies+=(
  "linear;\
  ${path_directory_gwas_concatenation}/testosterone_linear/male_unadjust_testosterone_imputation_log/gwas_concatenation.txt.gz;\
  ${path_directory_gwas_format}/ukbiobank_tcw_2022-05-04/testosterone_imputation_log_male_unadjust.txt.gz"
)

# within driver for-loop, call a script to call the appropriate format script...

# Organize information in format for LDSC.
for study_details in "${studies[@]}"; do
  # Read information.
  IFS=";" read -r -a array <<< "${study_details}"
  type_regression="${array[0]}"
  path_file_gwas_source="${array[1]}"
  path_file_gwas_product="${array[2]}"
  # Organize specific paths and parameters.
  report="true" # "true" or "false"
  if [[ "$type_regression" == "linear" ]]; then
    path_script_format_gwas="${path_script_format_gwas_linear}"
  elif [[ "$type_regression" == "logistic" ]]; then
    path_script_format_gwas="${path_script_format_gwas_logistic}"
  else
    echo "invalid specification of regression type"
  fi
  # Format adjustment.
  /usr/bin/bash "${path_script_format_gwas}" \
  $path_file_gwas_source \
  $path_file_gwas_product \
  $report
done
