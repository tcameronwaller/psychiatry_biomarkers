#!/bin/bash

###########################################################################
###########################################################################
###########################################################################
# ...
###########################################################################
###########################################################################
###########################################################################

# batch jobs: 132
# batch submission: TCW; at 13:17:20 on 10 August 2022

################################################################################
# Organize paths.
# Read private, local file paths.
cd ~/paths
path_gwas_summaries=$(<"./gwas_summaries_waller_metabolism.txt")
path_process=$(<"./process_psychiatric_metabolism.txt")
path_dock="$path_process/dock"
path_directory_gwas_source="${path_gwas_summaries}/tcameronwaller/ukbiobank_sex_steroid_hormones_proteins_2022-07-14"
path_directory_gwas_product="${path_dock}/gwas_format_team"

path_file_batch_instances="${path_directory_gwas_product}/batch_instances.txt"

# Scripts.
path_script_run_organize_format_gwas="$path_process/psychiatric_metabolism/scripts/record/2022-09-26/gwas_organization_team/2_run_batch_organize_format_gwas.sh"
path_promiscuity_scripts="${path_process}/promiscuity/scripts"
path_script_format_gwas_linear="${path_promiscuity_scripts}/gwas_process/format_gwas_team/format_gwas_team_plink_linear.sh"
path_script_format_gwas_logistic="${path_promiscuity_scripts}/gwas_process/format_gwas_team/format_gwas_team_plink_logistic.sh"

# Initialize directories.
rm -r $path_directory_gwas_product
mkdir -p $path_directory_gwas_product
rm $path_file_batch_instances

###########################################################################
# Execute procedure.

# Organize multi-dimensional array of information about studies.
studies=()
# [regression type] ; [full path to source file of GWAS summary statistics] ; [full path to product file of GWAS summary statistics]


##########

# Estradiol.
# oestradiol_imputation_log
# oestradiol_detection
studies+=(
  "linear;\
  ${path_directory_gwas_source}/oestradiol_linear/female_premenopause_joint_1_oestradiol_imputation_log/gwas.txt.gz;\
  ${path_directory_gwas_product}/tcw_ukb_estradiol_female_premenopause_adjust.txt.gz"
)
studies+=(
  "linear;\
  ${path_directory_gwas_source}/oestradiol_linear/female_premenopause_unadjust_oestradiol_imputation_log/gwas.txt.gz;\
  ${path_directory_gwas_product}/tcw_ukb_estradiol_female_premenopause_unadjust.txt.gz"
)
studies+=(
  "linear;\
  ${path_directory_gwas_source}/oestradiol_linear/female_perimenopause_joint_1_oestradiol_imputation_log/gwas.txt.gz;\
  ${path_directory_gwas_product}/tcw_ukb_estradiol_female_perimenopause_adjust.txt.gz"
)
studies+=(
  "linear;\
  ${path_directory_gwas_source}/oestradiol_linear/female_perimenopause_unadjust_oestradiol_imputation_log/gwas.txt.gz;\
  ${path_directory_gwas_product}/tcw_ukb_estradiol_female_perimenopause_unadjust.txt.gz"
)
studies+=(
  "logistic;\
  ${path_directory_gwas_source}/oestradiol_logistic/female_joint_1_oestradiol_detection/gwas.txt.gz;\
  ${path_directory_gwas_product}/tcw_ukb_estradiol_detection_female_adjust.txt.gz"
)
studies+=(
  "logistic;\
  ${path_directory_gwas_source}/oestradiol_logistic/female_unadjust_oestradiol_detection/gwas.txt.gz;\
  ${path_directory_gwas_product}/tcw_ukb_estradiol_detection_female_unadjust.txt.gz"
)
studies+=(
  "logistic;\
  ${path_directory_gwas_source}/oestradiol_logistic/female_premenopause_joint_1_oestradiol_detection/gwas.txt.gz;\
  ${path_directory_gwas_product}/tcw_ukb_estradiol_detection_female_premenopause_adjust.txt.gz"
)
studies+=(
  "logistic;\
  ${path_directory_gwas_source}/oestradiol_logistic/female_premenopause_unadjust_oestradiol_detection/gwas.txt.gz;\
  ${path_directory_gwas_product}/tcw_ukb_estradiol_detection_female_premenopause_unadjust.txt.gz"
)
studies+=(
  "logistic;\
  ${path_directory_gwas_source}/oestradiol_logistic/female_perimenopause_joint_1_oestradiol_detection/gwas.txt.gz;\
  ${path_directory_gwas_product}/tcw_ukb_estradiol_detection_female_perimenopause_adjust.txt.gz"
)
studies+=(
  "logistic;\
  ${path_directory_gwas_source}/oestradiol_logistic/female_perimenopause_unadjust_oestradiol_detection/gwas.txt.gz;\
  ${path_directory_gwas_product}/tcw_ukb_estradiol_detection_female_perimenopause_unadjust.txt.gz"
)
studies+=(
  "logistic;\
  ${path_directory_gwas_source}/oestradiol_logistic/female_postmenopause_joint_1_oestradiol_detection/gwas.txt.gz;\
  ${path_directory_gwas_product}/tcw_ukb_estradiol_detection_female_postmenopause_adjust.txt.gz"
)
studies+=(
  "logistic;\
  ${path_directory_gwas_source}/oestradiol_logistic/female_postmenopause_unadjust_oestradiol_detection/gwas.txt.gz;\
  ${path_directory_gwas_product}/tcw_ukb_estradiol_detection_female_postmenopause_unadjust.txt.gz"
)
studies+=(
  "logistic;\
  ${path_directory_gwas_source}/oestradiol_logistic/male_joint_1_oestradiol_detection/gwas.txt.gz;\
  ${path_directory_gwas_product}/tcw_ukb_estradiol_detection_male_adjust.txt.gz"
)
studies+=(
  "logistic;\
  ${path_directory_gwas_source}/oestradiol_logistic/male_unadjust_oestradiol_detection/gwas.txt.gz;\
  ${path_directory_gwas_product}/tcw_ukb_estradiol_detection_male_unadjust.txt.gz"
)
studies+=(
  "logistic;\
  ${path_directory_gwas_source}/oestradiol_logistic/male_age_low_joint_1_oestradiol_detection/gwas.txt.gz;\
  ${path_directory_gwas_product}/tcw_ukb_estradiol_detection_male_age_low_adjust.txt.gz"
)
studies+=(
  "logistic;\
  ${path_directory_gwas_source}/oestradiol_logistic/male_age_low_unadjust_oestradiol_detection/gwas.txt.gz;\
  ${path_directory_gwas_product}/tcw_ukb_estradiol_detection_male_age_low_unadjust.txt.gz"
)
studies+=(
  "logistic;\
  ${path_directory_gwas_source}/oestradiol_logistic/male_age_middle_joint_1_oestradiol_detection/gwas.txt.gz;\
  ${path_directory_gwas_product}/tcw_ukb_estradiol_detection_male_age_middle_adjust.txt.gz"
)
studies+=(
  "logistic;\
  ${path_directory_gwas_source}/oestradiol_logistic/male_age_middle_unadjust_oestradiol_detection/gwas.txt.gz;\
  ${path_directory_gwas_product}/tcw_ukb_estradiol_detection_male_age_middle_unadjust.txt.gz"
)
studies+=(
  "logistic;\
  ${path_directory_gwas_source}/oestradiol_logistic/male_age_high_joint_1_oestradiol_detection/gwas.txt.gz;\
  ${path_directory_gwas_product}/tcw_ukb_estradiol_detection_male_age_high_adjust.txt.gz"
)
studies+=(
  "logistic;\
  ${path_directory_gwas_source}/oestradiol_logistic/male_age_high_unadjust_oestradiol_detection/gwas.txt.gz;\
  ${path_directory_gwas_product}/tcw_ukb_estradiol_detection_male_age_high_unadjust.txt.gz"
)

# Estradiol, Bioavailable.
# oestradiol_bioavailable_imputation_log
studies+=(
  "linear;\
  ${path_directory_gwas_source}/oestradiol_bioavailable_linear/female_joint_1_oestradiol_bioavailable_imputation_log/gwas.txt.gz;\
  ${path_directory_gwas_product}/tcw_ukb_estradiol_bioavailable_female_adjust.txt.gz"
)
studies+=(
  "linear;\
  ${path_directory_gwas_source}/oestradiol_bioavailable_linear/female_unadjust_oestradiol_bioavailable_imputation_log/gwas.txt.gz;\
  ${path_directory_gwas_product}/tcw_ukb_estradiol_bioavailable_female_unadjust.txt.gz"
)
studies+=(
  "linear;\
  ${path_directory_gwas_source}/oestradiol_bioavailable_linear/female_premenopause_joint_1_oestradiol_bioavailable_imputation_log/gwas.txt.gz;\
  ${path_directory_gwas_product}/tcw_ukb_estradiol_bioavailable_female_premenopause_adjust.txt.gz"
)
studies+=(
  "linear;\
  ${path_directory_gwas_source}/oestradiol_bioavailable_linear/female_premenopause_unadjust_oestradiol_bioavailable_imputation_log/gwas.txt.gz;\
  ${path_directory_gwas_product}/tcw_ukb_estradiol_bioavailable_female_premenopause_unadjust.txt.gz"
)
studies+=(
  "linear;\
  ${path_directory_gwas_source}/oestradiol_bioavailable_linear/female_perimenopause_joint_1_oestradiol_bioavailable_imputation_log/gwas.txt.gz;\
  ${path_directory_gwas_product}/tcw_ukb_estradiol_bioavailable_female_perimenopause_adjust.txt.gz"
)
studies+=(
  "linear;\
  ${path_directory_gwas_source}/oestradiol_bioavailable_linear/female_perimenopause_unadjust_oestradiol_bioavailable_imputation_log/gwas.txt.gz;\
  ${path_directory_gwas_product}/tcw_ukb_estradiol_bioavailable_female_perimenopause_unadjust.txt.gz"
)
studies+=(
  "linear;\
  ${path_directory_gwas_source}/oestradiol_bioavailable_linear/female_postmenopause_joint_1_oestradiol_bioavailable_imputation_log/gwas.txt.gz;\
  ${path_directory_gwas_product}/tcw_ukb_estradiol_bioavailable_female_postmenopause_adjust.txt.gz"
)
studies+=(
  "linear;\
  ${path_directory_gwas_source}/oestradiol_bioavailable_linear/female_postmenopause_unadjust_oestradiol_bioavailable_imputation_log/gwas.txt.gz;\
  ${path_directory_gwas_product}/tcw_ukb_estradiol_bioavailable_female_postmenopause_unadjust.txt.gz"
)
studies+=(
  "linear;\
  ${path_directory_gwas_source}/oestradiol_bioavailable_linear/male_joint_1_oestradiol_bioavailable_imputation_log/gwas.txt.gz;\
  ${path_directory_gwas_product}/tcw_ukb_estradiol_bioavailable_male_adjust.txt.gz"
)
studies+=(
  "linear;\
  ${path_directory_gwas_source}/oestradiol_bioavailable_linear/male_unadjust_oestradiol_bioavailable_imputation_log/gwas.txt.gz;\
  ${path_directory_gwas_product}/tcw_ukb_estradiol_bioavailable_male_unadjust.txt.gz"
)
studies+=(
  "linear;\
  ${path_directory_gwas_source}/oestradiol_bioavailable_linear/male_age_low_joint_1_oestradiol_bioavailable_imputation_log/gwas.txt.gz;\
  ${path_directory_gwas_product}/tcw_ukb_estradiol_bioavailable_male_age_low_adjust.txt.gz"
)
studies+=(
  "linear;\
  ${path_directory_gwas_source}/oestradiol_bioavailable_linear/male_age_low_unadjust_oestradiol_bioavailable_imputation_log/gwas.txt.gz;\
  ${path_directory_gwas_product}/tcw_ukb_estradiol_bioavailable_male_age_low_unadjust.txt.gz"
)
studies+=(
  "linear;\
  ${path_directory_gwas_source}/oestradiol_bioavailable_linear/male_age_middle_joint_1_oestradiol_bioavailable_imputation_log/gwas.txt.gz;\
  ${path_directory_gwas_product}/tcw_ukb_estradiol_bioavailable_male_age_middle_adjust.txt.gz"
)
studies+=(
  "linear;\
  ${path_directory_gwas_source}/oestradiol_bioavailable_linear/male_age_middle_unadjust_oestradiol_bioavailable_imputation_log/gwas.txt.gz;\
  ${path_directory_gwas_product}/tcw_ukb_estradiol_bioavailable_male_age_middle_unadjust.txt.gz"
)
studies+=(
  "linear;\
  ${path_directory_gwas_source}/oestradiol_bioavailable_linear/male_age_high_joint_1_oestradiol_bioavailable_imputation_log/gwas.txt.gz;\
  ${path_directory_gwas_product}/tcw_ukb_estradiol_bioavailable_male_age_high_adjust.txt.gz"
)
studies+=(
  "linear;\
  ${path_directory_gwas_source}/oestradiol_bioavailable_linear/male_age_high_unadjust_oestradiol_bioavailable_imputation_log/gwas.txt.gz;\
  ${path_directory_gwas_product}/tcw_ukb_estradiol_bioavailable_male_age_high_unadjust.txt.gz"
)


# Estradiol, Free.
# oestradiol_free_imputation_log
studies+=(
  "linear;\
  ${path_directory_gwas_source}/oestradiol_free_linear/female_joint_1_oestradiol_free_imputation_log/gwas.txt.gz;\
  ${path_directory_gwas_product}/tcw_ukb_estradiol_free_female_adjust.txt.gz"
)
studies+=(
  "linear;\
  ${path_directory_gwas_source}/oestradiol_free_linear/female_unadjust_oestradiol_free_imputation_log/gwas.txt.gz;\
  ${path_directory_gwas_product}/tcw_ukb_estradiol_free_female_unadjust.txt.gz"
)
studies+=(
  "linear;\
  ${path_directory_gwas_source}/oestradiol_free_linear/female_premenopause_joint_1_oestradiol_free_imputation_log/gwas.txt.gz;\
  ${path_directory_gwas_product}/tcw_ukb_estradiol_free_female_premenopause_adjust.txt.gz"
)
studies+=(
  "linear;\
  ${path_directory_gwas_source}/oestradiol_free_linear/female_premenopause_unadjust_oestradiol_free_imputation_log/gwas.txt.gz;\
  ${path_directory_gwas_product}/tcw_ukb_estradiol_free_female_premenopause_unadjust.txt.gz"
)
studies+=(
  "linear;\
  ${path_directory_gwas_source}/oestradiol_free_linear/female_perimenopause_joint_1_oestradiol_free_imputation_log/gwas.txt.gz;\
  ${path_directory_gwas_product}/tcw_ukb_estradiol_free_female_perimenopause_adjust.txt.gz"
)
studies+=(
  "linear;\
  ${path_directory_gwas_source}/oestradiol_free_linear/female_perimenopause_unadjust_oestradiol_free_imputation_log/gwas.txt.gz;\
  ${path_directory_gwas_product}/tcw_ukb_estradiol_free_female_perimenopause_unadjust.txt.gz"
)
studies+=(
  "linear;\
  ${path_directory_gwas_source}/oestradiol_free_linear/female_postmenopause_joint_1_oestradiol_free_imputation_log/gwas.txt.gz;\
  ${path_directory_gwas_product}/tcw_ukb_estradiol_free_female_postmenopause_adjust.txt.gz"
)
studies+=(
  "linear;\
  ${path_directory_gwas_source}/oestradiol_free_linear/female_postmenopause_unadjust_oestradiol_free_imputation_log/gwas.txt.gz;\
  ${path_directory_gwas_product}/tcw_ukb_estradiol_free_female_postmenopause_unadjust.txt.gz"
)
studies+=(
  "linear;\
  ${path_directory_gwas_source}/oestradiol_free_linear/male_joint_1_oestradiol_free_imputation_log/gwas.txt.gz;\
  ${path_directory_gwas_product}/tcw_ukb_estradiol_free_male_adjust.txt.gz"
)
studies+=(
  "linear;\
  ${path_directory_gwas_source}/oestradiol_free_linear/male_unadjust_oestradiol_free_imputation_log/gwas.txt.gz;\
  ${path_directory_gwas_product}/tcw_ukb_estradiol_free_male_unadjust.txt.gz"
)
studies+=(
  "linear;\
  ${path_directory_gwas_source}/oestradiol_free_linear/male_age_low_joint_1_oestradiol_free_imputation_log/gwas.txt.gz;\
  ${path_directory_gwas_product}/tcw_ukb_estradiol_free_male_age_low_adjust.txt.gz"
)
studies+=(
  "linear;\
  ${path_directory_gwas_source}/oestradiol_free_linear/male_age_low_unadjust_oestradiol_free_imputation_log/gwas.txt.gz;\
  ${path_directory_gwas_product}/tcw_ukb_estradiol_free_male_age_low_unadjust.txt.gz"
)
studies+=(
  "linear;\
  ${path_directory_gwas_source}/oestradiol_free_linear/male_age_middle_joint_1_oestradiol_free_imputation_log/gwas.txt.gz;\
  ${path_directory_gwas_product}/tcw_ukb_estradiol_free_male_age_middle_adjust.txt.gz"
)
studies+=(
  "linear;\
  ${path_directory_gwas_source}/oestradiol_free_linear/male_age_middle_unadjust_oestradiol_free_imputation_log/gwas.txt.gz;\
  ${path_directory_gwas_product}/tcw_ukb_estradiol_free_male_age_middle_unadjust.txt.gz"
)
studies+=(
  "linear;\
  ${path_directory_gwas_source}/oestradiol_free_linear/male_age_high_joint_1_oestradiol_free_imputation_log/gwas.txt.gz;\
  ${path_directory_gwas_product}/tcw_ukb_estradiol_free_male_age_high_adjust.txt.gz"
)
studies+=(
  "linear;\
  ${path_directory_gwas_source}/oestradiol_free_linear/male_age_high_unadjust_oestradiol_free_imputation_log/gwas.txt.gz;\
  ${path_directory_gwas_product}/tcw_ukb_estradiol_free_male_age_high_unadjust.txt.gz"
)


##########

# Testosterone.
# testosterone_imputation_log
studies+=(
  "linear;\
  ${path_directory_gwas_source}/testosterone_linear/female_joint_1_testosterone_imputation_log/gwas.txt.gz;\
  ${path_directory_gwas_product}/tcw_ukb_testosterone_female_adjust.txt.gz"
)
studies+=(
  "linear;\
  ${path_directory_gwas_source}/testosterone_linear/female_unadjust_testosterone_imputation_log/gwas.txt.gz;\
  ${path_directory_gwas_product}/tcw_ukb_testosterone_female_unadjust.txt.gz"
)
studies+=(
  "linear;\
  ${path_directory_gwas_source}/testosterone_linear/female_premenopause_joint_1_testosterone_imputation_log/gwas.txt.gz;\
  ${path_directory_gwas_product}/tcw_ukb_testosterone_female_premenopause_adjust.txt.gz"
)
studies+=(
  "linear;\
  ${path_directory_gwas_source}/testosterone_linear/female_premenopause_unadjust_testosterone_imputation_log/gwas.txt.gz;\
  ${path_directory_gwas_product}/tcw_ukb_testosterone_female_premenopause_unadjust.txt.gz"
)
studies+=(
  "linear;\
  ${path_directory_gwas_source}/testosterone_linear/female_perimenopause_joint_1_testosterone_imputation_log/gwas.txt.gz;\
  ${path_directory_gwas_product}/tcw_ukb_testosterone_female_perimenopause_adjust.txt.gz"
)
studies+=(
  "linear;\
  ${path_directory_gwas_source}/testosterone_linear/female_perimenopause_unadjust_testosterone_imputation_log/gwas.txt.gz;\
  ${path_directory_gwas_product}/tcw_ukb_testosterone_female_perimenopause_unadjust.txt.gz"
)
studies+=(
  "linear;\
  ${path_directory_gwas_source}/testosterone_linear/female_postmenopause_joint_1_testosterone_imputation_log/gwas.txt.gz;\
  ${path_directory_gwas_product}/tcw_ukb_testosterone_female_postmenopause_adjust.txt.gz"
)
studies+=(
  "linear;\
  ${path_directory_gwas_source}/testosterone_linear/female_postmenopause_unadjust_testosterone_imputation_log/gwas.txt.gz;\
  ${path_directory_gwas_product}/tcw_ukb_testosterone_female_postmenopause_unadjust.txt.gz"
)
studies+=(
  "linear;\
  ${path_directory_gwas_source}/testosterone_linear/male_joint_1_testosterone_imputation_log/gwas.txt.gz;\
  ${path_directory_gwas_product}/tcw_ukb_testosterone_male_adjust.txt.gz"
)
studies+=(
  "linear;\
  ${path_directory_gwas_source}/testosterone_linear/male_unadjust_testosterone_imputation_log/gwas.txt.gz;\
  ${path_directory_gwas_product}/tcw_ukb_testosterone_male_unadjust.txt.gz"
)
studies+=(
  "linear;\
  ${path_directory_gwas_source}/testosterone_linear/male_age_low_joint_1_testosterone_imputation_log/gwas.txt.gz;\
  ${path_directory_gwas_product}/tcw_ukb_testosterone_male_age_low_adjust.txt.gz"
)
studies+=(
  "linear;\
  ${path_directory_gwas_source}/testosterone_linear/male_age_low_unadjust_testosterone_imputation_log/gwas.txt.gz;\
  ${path_directory_gwas_product}/tcw_ukb_testosterone_male_age_low_unadjust.txt.gz"
)
studies+=(
  "linear;\
  ${path_directory_gwas_source}/testosterone_linear/male_age_middle_joint_1_testosterone_imputation_log/gwas.txt.gz;\
  ${path_directory_gwas_product}/tcw_ukb_testosterone_male_age_middle_adjust.txt.gz"
)
studies+=(
  "linear;\
  ${path_directory_gwas_source}/testosterone_linear/male_age_middle_unadjust_testosterone_imputation_log/gwas.txt.gz;\
  ${path_directory_gwas_product}/tcw_ukb_testosterone_male_age_middle_unadjust.txt.gz"
)
studies+=(
  "linear;\
  ${path_directory_gwas_source}/testosterone_linear/male_age_high_joint_1_testosterone_imputation_log/gwas.txt.gz;\
  ${path_directory_gwas_product}/tcw_ukb_testosterone_male_age_high_adjust.txt.gz"
)
studies+=(
  "linear;\
  ${path_directory_gwas_source}/testosterone_linear/male_age_high_unadjust_testosterone_imputation_log/gwas.txt.gz;\
  ${path_directory_gwas_product}/tcw_ukb_testosterone_male_age_high_unadjust.txt.gz"
)

# Testosterone, Bioavailable.
# testosterone_bioavailable_imputation_log
studies+=(
  "linear;\
  ${path_directory_gwas_source}/testosterone_bioavailable_linear/female_joint_1_testosterone_bioavailable_imputation_log/gwas.txt.gz;\
  ${path_directory_gwas_product}/tcw_ukb_testosterone_bioavailable_female_adjust.txt.gz"
)
studies+=(
  "linear;\
  ${path_directory_gwas_source}/testosterone_bioavailable_linear/female_unadjust_testosterone_bioavailable_imputation_log/gwas.txt.gz;\
  ${path_directory_gwas_product}/tcw_ukb_testosterone_bioavailable_female_unadjust.txt.gz"
)
studies+=(
  "linear;\
  ${path_directory_gwas_source}/testosterone_bioavailable_linear/female_premenopause_joint_1_testosterone_bioavailable_imputation_log/gwas.txt.gz;\
  ${path_directory_gwas_product}/tcw_ukb_testosterone_bioavailable_female_premenopause_adjust.txt.gz"
)
studies+=(
  "linear;\
  ${path_directory_gwas_source}/testosterone_bioavailable_linear/female_premenopause_unadjust_testosterone_bioavailable_imputation_log/gwas.txt.gz;\
  ${path_directory_gwas_product}/tcw_ukb_testosterone_bioavailable_female_premenopause_unadjust.txt.gz"
)
studies+=(
  "linear;\
  ${path_directory_gwas_source}/testosterone_bioavailable_linear/female_perimenopause_joint_1_testosterone_bioavailable_imputation_log/gwas.txt.gz;\
  ${path_directory_gwas_product}/tcw_ukb_testosterone_bioavailable_female_perimenopause_adjust.txt.gz"
)
studies+=(
  "linear;\
  ${path_directory_gwas_source}/testosterone_bioavailable_linear/female_perimenopause_unadjust_testosterone_bioavailable_imputation_log/gwas.txt.gz;\
  ${path_directory_gwas_product}/tcw_ukb_testosterone_bioavailable_female_perimenopause_unadjust.txt.gz"
)
studies+=(
  "linear;\
  ${path_directory_gwas_source}/testosterone_bioavailable_linear/female_postmenopause_joint_1_testosterone_bioavailable_imputation_log/gwas.txt.gz;\
  ${path_directory_gwas_product}/tcw_ukb_testosterone_bioavailable_female_postmenopause_adjust.txt.gz"
)
studies+=(
  "linear;\
  ${path_directory_gwas_source}/testosterone_bioavailable_linear/female_postmenopause_unadjust_testosterone_bioavailable_imputation_log/gwas.txt.gz;\
  ${path_directory_gwas_product}/tcw_ukb_testosterone_bioavailable_female_postmenopause_unadjust.txt.gz"
)
studies+=(
  "linear;\
  ${path_directory_gwas_source}/testosterone_bioavailable_linear/male_joint_1_testosterone_bioavailable_imputation_log/gwas.txt.gz;\
  ${path_directory_gwas_product}/tcw_ukb_testosterone_bioavailable_male_adjust.txt.gz"
)
studies+=(
  "linear;\
  ${path_directory_gwas_source}/testosterone_bioavailable_linear/male_unadjust_testosterone_bioavailable_imputation_log/gwas.txt.gz;\
  ${path_directory_gwas_product}/tcw_ukb_testosterone_bioavailable_male_unadjust.txt.gz"
)
studies+=(
  "linear;\
  ${path_directory_gwas_source}/testosterone_bioavailable_linear/male_age_low_joint_1_testosterone_bioavailable_imputation_log/gwas.txt.gz;\
  ${path_directory_gwas_product}/tcw_ukb_testosterone_bioavailable_male_age_low_adjust.txt.gz"
)
studies+=(
  "linear;\
  ${path_directory_gwas_source}/testosterone_bioavailable_linear/male_age_low_unadjust_testosterone_bioavailable_imputation_log/gwas.txt.gz;\
  ${path_directory_gwas_product}/tcw_ukb_testosterone_bioavailable_male_age_low_unadjust.txt.gz"
)
studies+=(
  "linear;\
  ${path_directory_gwas_source}/testosterone_bioavailable_linear/male_age_middle_joint_1_testosterone_bioavailable_imputation_log/gwas.txt.gz;\
  ${path_directory_gwas_product}/tcw_ukb_testosterone_bioavailable_male_age_middle_adjust.txt.gz"
)
studies+=(
  "linear;\
  ${path_directory_gwas_source}/testosterone_bioavailable_linear/male_age_middle_unadjust_testosterone_bioavailable_imputation_log/gwas.txt.gz;\
  ${path_directory_gwas_product}/tcw_ukb_testosterone_bioavailable_male_age_middle_unadjust.txt.gz"
)
studies+=(
  "linear;\
  ${path_directory_gwas_source}/testosterone_bioavailable_linear/male_age_high_joint_1_testosterone_bioavailable_imputation_log/gwas.txt.gz;\
  ${path_directory_gwas_product}/tcw_ukb_testosterone_bioavailable_male_age_high_adjust.txt.gz"
)
studies+=(
  "linear;\
  ${path_directory_gwas_source}/testosterone_bioavailable_linear/male_age_high_unadjust_testosterone_bioavailable_imputation_log/gwas.txt.gz;\
  ${path_directory_gwas_product}/tcw_ukb_testosterone_bioavailable_male_age_high_unadjust.txt.gz"
)

# Testosterone, Free.
# testosterone_free_imputation_log
studies+=(
  "linear;\
  ${path_directory_gwas_source}/testosterone_free_linear/female_joint_1_testosterone_free_imputation_log/gwas.txt.gz;\
  ${path_directory_gwas_product}/tcw_ukb_testosterone_free_female_adjust.txt.gz"
)
studies+=(
  "linear;\
  ${path_directory_gwas_source}/testosterone_free_linear/female_unadjust_testosterone_free_imputation_log/gwas.txt.gz;\
  ${path_directory_gwas_product}/tcw_ukb_testosterone_free_female_unadjust.txt.gz"
)
studies+=(
  "linear;\
  ${path_directory_gwas_source}/testosterone_free_linear/female_premenopause_joint_1_testosterone_free_imputation_log/gwas.txt.gz;\
  ${path_directory_gwas_product}/tcw_ukb_testosterone_free_female_premenopause_adjust.txt.gz"
)
studies+=(
  "linear;\
  ${path_directory_gwas_source}/testosterone_free_linear/female_premenopause_unadjust_testosterone_free_imputation_log/gwas.txt.gz;\
  ${path_directory_gwas_product}/tcw_ukb_testosterone_free_female_premenopause_unadjust.txt.gz"
)
studies+=(
  "linear;\
  ${path_directory_gwas_source}/testosterone_free_linear/female_perimenopause_joint_1_testosterone_free_imputation_log/gwas.txt.gz;\
  ${path_directory_gwas_product}/tcw_ukb_testosterone_free_female_perimenopause_adjust.txt.gz"
)
studies+=(
  "linear;\
  ${path_directory_gwas_source}/testosterone_free_linear/female_perimenopause_unadjust_testosterone_free_imputation_log/gwas.txt.gz;\
  ${path_directory_gwas_product}/tcw_ukb_testosterone_free_female_perimenopause_unadjust.txt.gz"
)
studies+=(
  "linear;\
  ${path_directory_gwas_source}/testosterone_free_linear/female_postmenopause_joint_1_testosterone_free_imputation_log/gwas.txt.gz;\
  ${path_directory_gwas_product}/tcw_ukb_testosterone_free_female_postmenopause_adjust.txt.gz"
)
studies+=(
  "linear;\
  ${path_directory_gwas_source}/testosterone_free_linear/female_postmenopause_unadjust_testosterone_free_imputation_log/gwas.txt.gz;\
  ${path_directory_gwas_product}/tcw_ukb_testosterone_free_female_postmenopause_unadjust.txt.gz"
)
studies+=(
  "linear;\
  ${path_directory_gwas_source}/testosterone_free_linear/male_joint_1_testosterone_free_imputation_log/gwas.txt.gz;\
  ${path_directory_gwas_product}/tcw_ukb_testosterone_free_male_adjust.txt.gz"
)
studies+=(
  "linear;\
  ${path_directory_gwas_source}/testosterone_free_linear/male_unadjust_testosterone_free_imputation_log/gwas.txt.gz;\
  ${path_directory_gwas_product}/tcw_ukb_testosterone_free_male_unadjust.txt.gz"
)
studies+=(
  "linear;\
  ${path_directory_gwas_source}/testosterone_free_linear/male_age_low_joint_1_testosterone_free_imputation_log/gwas.txt.gz;\
  ${path_directory_gwas_product}/tcw_ukb_testosterone_free_male_age_low_adjust.txt.gz"
)
studies+=(
  "linear;\
  ${path_directory_gwas_source}/testosterone_free_linear/male_age_low_unadjust_testosterone_free_imputation_log/gwas.txt.gz;\
  ${path_directory_gwas_product}/tcw_ukb_testosterone_free_male_age_low_unadjust.txt.gz"
)
studies+=(
  "linear;\
  ${path_directory_gwas_source}/testosterone_free_linear/male_age_middle_joint_1_testosterone_free_imputation_log/gwas.txt.gz;\
  ${path_directory_gwas_product}/tcw_ukb_testosterone_free_male_age_middle_adjust.txt.gz"
)
studies+=(
  "linear;\
  ${path_directory_gwas_source}/testosterone_free_linear/male_age_middle_unadjust_testosterone_free_imputation_log/gwas.txt.gz;\
  ${path_directory_gwas_product}/tcw_ukb_testosterone_free_male_age_middle_unadjust.txt.gz"
)
studies+=(
  "linear;\
  ${path_directory_gwas_source}/testosterone_free_linear/male_age_high_joint_1_testosterone_free_imputation_log/gwas.txt.gz;\
  ${path_directory_gwas_product}/tcw_ukb_testosterone_free_male_age_high_adjust.txt.gz"
)
studies+=(
  "linear;\
  ${path_directory_gwas_source}/testosterone_free_linear/male_age_high_unadjust_testosterone_free_imputation_log/gwas.txt.gz;\
  ${path_directory_gwas_product}/tcw_ukb_testosterone_free_male_age_high_unadjust.txt.gz"
)

##########

# Sex Steroid Hormone Binding Globulin (SHBG).
# steroid_globulin_imputation_log
studies+=(
  "linear;\
  ${path_directory_gwas_source}/steroid_globulin_linear/female_joint_1_steroid_globulin_imputation_log/gwas.txt.gz;\
  ${path_directory_gwas_product}/tcw_ukb_shbg_female_adjust.txt.gz"
)
studies+=(
  "linear;\
  ${path_directory_gwas_source}/steroid_globulin_linear/female_unadjust_steroid_globulin_imputation_log/gwas.txt.gz;\
  ${path_directory_gwas_product}/tcw_ukb_shbg_female_unadjust.txt.gz"
)
studies+=(
  "linear;\
  ${path_directory_gwas_source}/steroid_globulin_linear/female_premenopause_joint_1_steroid_globulin_imputation_log/gwas.txt.gz;\
  ${path_directory_gwas_product}/tcw_ukb_shbg_female_premenopause_adjust.txt.gz"
)
studies+=(
  "linear;\
  ${path_directory_gwas_source}/steroid_globulin_linear/female_premenopause_unadjust_steroid_globulin_imputation_log/gwas.txt.gz;\
  ${path_directory_gwas_product}/tcw_ukb_shbg_female_premenopause_unadjust.txt.gz"
)
studies+=(
  "linear;\
  ${path_directory_gwas_source}/steroid_globulin_linear/female_perimenopause_joint_1_steroid_globulin_imputation_log/gwas.txt.gz;\
  ${path_directory_gwas_product}/tcw_ukb_shbg_female_perimenopause_adjust.txt.gz"
)
studies+=(
  "linear;\
  ${path_directory_gwas_source}/steroid_globulin_linear/female_perimenopause_unadjust_steroid_globulin_imputation_log/gwas.txt.gz;\
  ${path_directory_gwas_product}/tcw_ukb_shbg_female_perimenopause_unadjust.txt.gz"
)
studies+=(
  "linear;\
  ${path_directory_gwas_source}/steroid_globulin_linear/female_postmenopause_joint_1_steroid_globulin_imputation_log/gwas.txt.gz;\
  ${path_directory_gwas_product}/tcw_ukb_shbg_female_postmenopause_adjust.txt.gz"
)
studies+=(
  "linear;\
  ${path_directory_gwas_source}/steroid_globulin_linear/female_postmenopause_unadjust_steroid_globulin_imputation_log/gwas.txt.gz;\
  ${path_directory_gwas_product}/tcw_ukb_shbg_female_postmenopause_unadjust.txt.gz"
)
studies+=(
  "linear;\
  ${path_directory_gwas_source}/steroid_globulin_linear/male_joint_1_steroid_globulin_imputation_log/gwas.txt.gz;\
  ${path_directory_gwas_product}/tcw_ukb_shbg_male_adjust.txt.gz"
)
studies+=(
  "linear;\
  ${path_directory_gwas_source}/steroid_globulin_linear/male_unadjust_steroid_globulin_imputation_log/gwas.txt.gz;\
  ${path_directory_gwas_product}/tcw_ukb_shbg_male_unadjust.txt.gz"
)
studies+=(
  "linear;\
  ${path_directory_gwas_source}/steroid_globulin_linear/male_age_low_joint_1_steroid_globulin_imputation_log/gwas.txt.gz;\
  ${path_directory_gwas_product}/tcw_ukb_shbg_male_age_low_adjust.txt.gz"
)
studies+=(
  "linear;\
  ${path_directory_gwas_source}/steroid_globulin_linear/male_age_low_unadjust_steroid_globulin_imputation_log/gwas.txt.gz;\
  ${path_directory_gwas_product}/tcw_ukb_shbg_male_age_low_unadjust.txt.gz"
)
studies+=(
  "linear;\
  ${path_directory_gwas_source}/steroid_globulin_linear/male_age_middle_joint_1_steroid_globulin_imputation_log/gwas.txt.gz;\
  ${path_directory_gwas_product}/tcw_ukb_shbg_male_age_middle_adjust.txt.gz"
)
studies+=(
  "linear;\
  ${path_directory_gwas_source}/steroid_globulin_linear/male_age_middle_unadjust_steroid_globulin_imputation_log/gwas.txt.gz;\
  ${path_directory_gwas_product}/tcw_ukb_shbg_male_age_middle_unadjust.txt.gz"
)
studies+=(
  "linear;\
  ${path_directory_gwas_source}/steroid_globulin_linear/male_age_high_joint_1_steroid_globulin_imputation_log/gwas.txt.gz;\
  ${path_directory_gwas_product}/tcw_ukb_shbg_male_age_high_adjust.txt.gz"
)
studies+=(
  "linear;\
  ${path_directory_gwas_source}/steroid_globulin_linear/male_age_high_unadjust_steroid_globulin_imputation_log/gwas.txt.gz;\
  ${path_directory_gwas_product}/tcw_ukb_shbg_male_age_high_unadjust.txt.gz"
)

##########

# Albumin.
# albumin_imputation
studies+=(
  "linear;\
  ${path_directory_gwas_source}/albumin_linear/female_joint_1_albumin_imputation/gwas.txt.gz;\
  ${path_directory_gwas_product}/tcw_ukb_albumin_female_adjust.txt.gz"
)
studies+=(
  "linear;\
  ${path_directory_gwas_source}/albumin_linear/female_unadjust_albumin_imputation/gwas.txt.gz;\
  ${path_directory_gwas_product}/tcw_ukb_albumin_female_unadjust.txt.gz"
)
studies+=(
  "linear;\
  ${path_directory_gwas_source}/albumin_linear/female_premenopause_joint_1_albumin_imputation/gwas.txt.gz;\
  ${path_directory_gwas_product}/tcw_ukb_albumin_female_premenopause_adjust.txt.gz"
)
studies+=(
  "linear;\
  ${path_directory_gwas_source}/albumin_linear/female_premenopause_unadjust_albumin_imputation/gwas.txt.gz;\
  ${path_directory_gwas_product}/tcw_ukb_albumin_female_premenopause_unadjust.txt.gz"
)
studies+=(
  "linear;\
  ${path_directory_gwas_source}/albumin_linear/female_perimenopause_joint_1_albumin_imputation/gwas.txt.gz;\
  ${path_directory_gwas_product}/tcw_ukb_albumin_female_perimenopause_adjust.txt.gz"
)
studies+=(
  "linear;\
  ${path_directory_gwas_source}/albumin_linear/female_perimenopause_unadjust_albumin_imputation/gwas.txt.gz;\
  ${path_directory_gwas_product}/tcw_ukb_albumin_female_perimenopause_unadjust.txt.gz"
)
studies+=(
  "linear;\
  ${path_directory_gwas_source}/albumin_linear/female_postmenopause_joint_1_albumin_imputation/gwas.txt.gz;\
  ${path_directory_gwas_product}/tcw_ukb_albumin_female_postmenopause_adjust.txt.gz"
)
studies+=(
  "linear;\
  ${path_directory_gwas_source}/albumin_linear/female_postmenopause_unadjust_albumin_imputation/gwas.txt.gz;\
  ${path_directory_gwas_product}/tcw_ukb_albumin_female_postmenopause_unadjust.txt.gz"
)
studies+=(
  "linear;\
  ${path_directory_gwas_source}/albumin_linear/male_joint_1_albumin_imputation/gwas.txt.gz;\
  ${path_directory_gwas_product}/tcw_ukb_albumin_male_adjust.txt.gz"
)
studies+=(
  "linear;\
  ${path_directory_gwas_source}/albumin_linear/male_unadjust_albumin_imputation/gwas.txt.gz;\
  ${path_directory_gwas_product}/tcw_ukb_albumin_male_unadjust.txt.gz"
)
studies+=(
  "linear;\
  ${path_directory_gwas_source}/albumin_linear/male_age_low_joint_1_albumin_imputation/gwas.txt.gz;\
  ${path_directory_gwas_product}/tcw_ukb_albumin_male_age_low_adjust.txt.gz"
)
studies+=(
  "linear;\
  ${path_directory_gwas_source}/albumin_linear/male_age_low_unadjust_albumin_imputation/gwas.txt.gz;\
  ${path_directory_gwas_product}/tcw_ukb_albumin_male_age_low_unadjust.txt.gz"
)
studies+=(
  "linear;\
  ${path_directory_gwas_source}/albumin_linear/male_age_middle_joint_1_albumin_imputation/gwas.txt.gz;\
  ${path_directory_gwas_product}/tcw_ukb_albumin_male_age_middle_adjust.txt.gz"
)
studies+=(
  "linear;\
  ${path_directory_gwas_source}/albumin_linear/male_age_middle_unadjust_albumin_imputation/gwas.txt.gz;\
  ${path_directory_gwas_product}/tcw_ukb_albumin_male_age_middle_unadjust.txt.gz"
)
studies+=(
  "linear;\
  ${path_directory_gwas_source}/albumin_linear/male_age_high_joint_1_albumin_imputation/gwas.txt.gz;\
  ${path_directory_gwas_product}/tcw_ukb_albumin_male_age_high_adjust.txt.gz"
)
studies+=(
  "linear;\
  ${path_directory_gwas_source}/albumin_linear/male_age_high_unadjust_albumin_imputation/gwas.txt.gz;\
  ${path_directory_gwas_product}/tcw_ukb_albumin_male_age_high_unadjust.txt.gz"
)

################################################################################

# Organize information for batch instances.
for study_details in "${studies[@]}"; do
  # Separate fields from instance.
  # [regression type] ; [full path to source file of GWAS summary statistics] ; [full path to product file of GWAS summary statistics]
  IFS=";" read -r -a array <<< "${study_details}"
  type_regression="${array[0]}"
  path_file_gwas_source="${array[1]}"
  path_file_gwas_product="${array[2]}"
  # Define and append a new batch instance.
  instance="${type_regression};${path_file_gwas_source};${path_file_gwas_product}"
  echo $instance >> $path_file_batch_instances
done

################################################################################
# Report.
if [[ "$report" == "true" ]]; then
  echo "----------"
  echo "1_submit_batch_organize_format_gwas.sh"
  echo "----------"
fi

################################################################################
# Submit batch instances to cluster scheduler.

# Read batch instances.
readarray -t batch_instances < $path_file_batch_instances
batch_instances_count=${#batch_instances[@]}
echo "----------"
echo "count of batch instances: " $batch_instances_count
echo "first batch instance: " ${batch_instances[0]} # notice base-zero indexing
echo "last batch instance: " ${batch_instances[$batch_instances_count - 1]}

# Execute batch with grid scheduler.
if true; then
  # Submit array batch to Sun Grid Engine.
  # Array batch indices must start at one (not zero).
  qsub -t 1-${batch_instances_count}:1 -o \
  "${path_directory_gwas_product}/batch_out.txt" -e "${path_directory_gwas_product}/batch_error.txt" \
  "${path_script_run_organize_format_gwas}" \
  $path_file_batch_instances \
  $batch_instances_count \
  $path_script_format_gwas_linear \
  $path_script_format_gwas_logistic \
  $report
fi



#
