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
path_batch_instances="${path_directory_gwas_format}/batch_instances.txt"

# Scripts.
path_script_run_organize_format_gwas="$path_process/psychiatric_metabolism/scripts/record/2022-05-04/gwas_organization/2_run_batch_organize_format_gwas.sh"
path_promiscuity_scripts="${path_process}/promiscuity/scripts"
path_script_format_gwas_linear="${path_promiscuity_scripts}/gwas_process/format_gwas_team/format_gwas_team_plink_linear.sh"
path_script_format_gwas_logistic="${path_promiscuity_scripts}/gwas_process/format_gwas_team/format_gwas_team_plink_logistic.sh"

# Initialize directories.
rm -r $path_directory_gwas_format
mkdir -p $path_directory_gwas_format
rm $path_batch_instances

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

# Organize information for batch instances.
for study_details in "${studies[@]}"; do
  # Separate fields from instance.
  IFS=";" read -r -a array <<< "${study_details}"
  type_regression="${array[0]}"
  path_file_gwas_source="${array[1]}"
  path_file_gwas_product="${array[2]}"
  # Define and append a new batch instance.
  instance="${type_regression};${path_file_gwas_source};${path_file_gwas_product}"
  echo $instance >> $path_batch_instances
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
readarray -t batch_instances < $path_batch_instances
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
  "${path_directory_gwas_format}/batch_out.txt" -e "${path_directory_gwas_format}/batch_error.txt" \
  "${path_script_run_organize_format_gwas}" \
  $path_batch_instances \
  $batch_instances_count \
  $path_script_format_gwas_linear \
  $path_script_format_gwas_logistic \
  $report
fi



#
