#!/bin/bash

###########################################################################
###########################################################################
###########################################################################
# ...
###########################################################################
###########################################################################
###########################################################################


################################################################################
# General parameters.

#cohorts_models="albumin_linear_2"                      #  2 GWAS; __ April 2022
#cohorts_models="albumin_linear_1"                      # 18 GWAS; __ April 2022
#cohorts_models="steroid_globulin_linear_2"             #  2 GWAS; __ April 2022
#cohorts_models="steroid_globulin_linear_1"             # 18 GWAS; __ April 2022
#cohorts_models="oestradiol_logistic"                   # __ GWAS; incomplete
#cohorts_models="oestradiol_linear_1"                   # 18 GWAS; __ April 2022
#cohorts_models="oestradiol_linear_2"                   # 12 GWAS; __ April 2022
#cohorts_models="oestradiol_bioavailable_linear"        # 18 GWAS; __ April 2022
#cohorts_models="oestradiol_free_linear"                # 18 GWAS; __ April 2022
#cohorts_models="testosterone_logistic"                 # __ GWAS; incomplete
#cohorts_models="testosterone_linear"                   # 54 GWAS; __ April 2022
#cohorts_models="testosterone_bioavailable_linear"      # 18 GWAS; __ April 2022
#cohorts_models="testosterone_free_linear"              # 18 GWAS; __ April 2022

# Parameters.
regression_type="linear" # "linear" or "logistic"
response="coefficient" # "coefficient", "odds_ratio", or "z_score"; for linear GWAS, use "coefficient" unless "response_standard_scale" is "true", in which case "z_score"
response_standard_scale="false" # whether to convert reponse (effect, coefficient) to z-score standard scale ("true" or "false")

restore_target_study_directories="true" # whether to delete any previous directories for each study's format and munge GWAS ("true" or "false")

################################################################################
# Organize paths.
# Read private, local file paths.
cd ~/paths
path_process=$(<"./process_psychiatric_metabolism.txt")
path_dock="$path_process/dock"

path_gwas_concatenation_container="${path_dock}/gwas_concatenation/${cohorts_models}"
path_gwas_format_container="${path_dock}/gwas_ldsc_format/${cohorts_models}"
path_gwas_munge_container="${path_dock}/gwas_ldsc_munge/${cohorts_models}"
path_heritability_container="${path_dock}/heritability/${cohorts_models}"

path_scripts_record="$path_process/psychiatric_metabolism/scripts/record/2022-05-04/ldsc_heritability_correlation"
path_batch_instances="${path_gwas_munge_container}/batch_instances_format_munge.txt"

###########################################################################
# Define explicit inclusions and exclusions.
# Use inclusions to run procedure for a few specific cohort-hormone combinations that are missing from the set.
# Use exclusions to omit a few cohort-hormone combinations that are not complete yet.

#delimiter=" "
#IFS=${delimiter}
#exclusions=()
#exclusions+=("female_combination_unadjust_albumin_log")
#unset IFS

###########################################################################
# Execute procedure.

# Initialize directories and batch instances.
rm -r $path_gwas_format_container
mkdir -p $path_gwas_format_container
rm -r $path_gwas_munge_container
mkdir -p $path_gwas_munge_container
rm -r $path_heritability_container
mkdir -p $path_heritability_container
rm $path_batch_instances

# Iterate on directories for GWAS on cohorts and hormones.
name_gwas_concatenation_file="gwas_concatenation.txt.gz"
cd $path_gwas_concatenation_container
for path_directory in `find . -maxdepth 1 -mindepth 1 -type d -not -name .`; do
  if [ -d "$path_directory" ]; then
    # Current content item is a directory.
    # Extract directory's base name.
    study="$(basename -- $path_directory)"
    #echo $directory
    # Determine whether directory contains valid GWAS summary statistics.
    matches=$(find "${path_gwas_concatenation_container}/${study}" -name "$name_gwas_concatenation_file")
    match_file=${matches[0]}
    if [[ -n $matches && -f $match_file ]]; then
      instance="$study;${path_gwas_concatenation_container}/${study}/${name_gwas_concatenation_file}"
      echo $instance >> $path_batch_instances
    fi
  fi
done

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
  "${path_gwas_munge_container}/post_process_out.txt" -e "${path_gwas_munge_container}/post_process_error.txt" \
  "${path_scripts_record}/4_run_batch_jobs_gwas_format_munge_heritability.sh" \
  $path_batch_instances \
  $batch_instances_count \
  $regression_type \
  $response \
  $response_standard_scale \
  $path_gwas_format_container \
  $path_gwas_munge_container \
  $path_heritability_container \
  $path_scripts_record \
  $path_process \
  $restore_target_study_directories
fi
