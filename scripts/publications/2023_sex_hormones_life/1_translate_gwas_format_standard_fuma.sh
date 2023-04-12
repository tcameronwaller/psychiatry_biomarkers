#!/bin/bash

################################################################################
# Author: T. Cameron Waller
# Date, first execution: 12 April 2023
# Date, last execution: 12 April 2023
# Review: TCW; 12 April 2023
################################################################################
# Note



################################################################################
# Organize paths.

# Directories.
cd ~/paths
path_directory_summaries=$(<"./gwas_summaries_waller_metabolism.txt")
path_directory_process=$(<"./process_psychiatric_metabolism.txt")
path_directory_dock="${path_directory_process}/dock"
path_directory_source="${path_directory_summaries}/tcameronwaller/ukbiobank_sex_steroid_hormones_proteins_2022-07-14"
path_directory_product="${path_directory_dock}/gwas_tcw_hormones_fuma"

# Scripts.
path_directory_promiscuity_scripts="${path_directory_process}/promiscuity/scripts"
path_file_script_standard="${path_directory_promiscuity_scripts}/gwas_format/translate_format_raw_to_standard/translate_gwas_plink_linear.sh"
path_file_script_fuma="${path_directory_promiscuity_scripts}/fuma/translate_gwas_standard_to_fuma.sh"

# Initialize directories.
rm -r $path_directory_product
mkdir -p $path_directory_product

################################################################################
# Organize parameters.

report="true"

################################################################################
# Execute procedure.

# Organize multi-dimensional array of information about studies.
studies=()
# [regression type] ; [full path to source file of GWAS summary statistics] ; [full path to product file of GWAS summary statistics]

##########

# Testosterone.
# testosterone_imputation_log
studies+=(
  "linear;\
  ${path_directory_source}/testosterone_linear/female_premenopause_joint_1_testosterone_imputation_log/gwas.txt.gz;\
  ${path_directory_product}/tcw_ukb_testosterone_female_premenopause_adjust_standard.txt.gz;\
  ${path_directory_product}/tcw_ukb_testosterone_female_premenopause_adjust_fuma.txt.gz"
)
studies+=(
  "linear;\
  ${path_directory_source}/testosterone_linear/female_premenopause_unadjust_testosterone_imputation_log/gwas.txt.gz;\
  ${path_directory_product}/tcw_ukb_testosterone_female_premenopause_unadjust_standard.txt.gz;\
  ${path_directory_product}/tcw_ukb_testosterone_female_premenopause_unadjust_fuma.txt.gz"
)
studies+=(
  "linear;\
  ${path_directory_source}/testosterone_linear/female_postmenopause_joint_1_testosterone_imputation_log/gwas.txt.gz;\
  ${path_directory_product}/tcw_ukb_testosterone_female_postmenopause_adjust_standard.txt.gz;\
  ${path_directory_product}/tcw_ukb_testosterone_female_postmenopause_adjust_fuma.txt.gz"
)
studies+=(
  "linear;\
  ${path_directory_source}/testosterone_linear/female_postmenopause_unadjust_testosterone_imputation_log/gwas.txt.gz;\
  ${path_directory_product}/tcw_ukb_testosterone_female_postmenopause_unadjust_standard.txt.gz;\
  ${path_directory_product}/tcw_ukb_testosterone_female_postmenopause_unadjust_fuma.txt.gz"
)


################################################################################

# Organize information for batch instances.
for study_details in "${studies[@]}"; do
  # Separate fields from instance.
  # [regression type] ; [full path to source file of GWAS summary statistics] ; [full path to product file of GWAS summary statistics]
  IFS=";" read -r -a array <<< "${study_details}"
  type_regression="${array[0]}"
  path_file_source="${array[1]}"
  path_file_product_standard="${array[2]}"
  path_file_product_fuma="${array[3]}"

  # Translate format to team standard.
  /usr/bin/bash $path_file_script_standard \
  $path_file_source \
  $path_file_product_standard \
  "0" \
  "NA" \
  "0" \
  "NA" \
  "NA" \
  $report

  # Translate format to abbreviation for upload to FUMA.
  /usr/bin/bash $path_file_script_fuma \
  $path_file_product_standard \
  $path_file_product_fuma \
  $report

done

################################################################################
# Report.
if [[ "$report" == "true" ]]; then
  echo "----------"
  echo "Script complete:"
  echo $0
  echo "1_translate_gwas_format_standard_fuma.sh"
  echo "----------"
fi

#
