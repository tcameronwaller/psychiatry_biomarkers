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

cohorts_models="body_white_bipolar_strict"
#cohorts_models="body_white_bipolar_loose"

name_gwas_munge_file="gwas_munge.sumstats.gz"

################################################################################
# Organize paths.
# Read private, local file paths.
cd ~/paths
path_process=$(<"./process_psychiatric_metabolism.txt")
path_dock="$path_process/dock"

path_gwas_source_container="${path_dock}/gwas_process/${cohorts_models}"

###########################################################################
# Define main comparisons.

# Define array of primary studies.
primaries=()

primaries+=("30124842_yengo_2018;${path_dock}/gwas_process/30124842_yengo_2018/${name_gwas_munge_file}")

# Define array of secondary studies.
secondaries=()
# Iterate on directories for GWAS on cohorts and hormones.
cd $path_gwas_source_container
for path_directory in `find . -maxdepth 1 -mindepth 1 -type d -not -name .`; do
  if [ -d "$path_directory" ]; then
    # Current content item is a directory.
    # Extract directory's base name.
    study="$(basename -- $path_directory)"
    #echo $directory
    # Determine whether directory contains valid GWAS summary statistics.
    matches=$(find "${path_gwas_source_container}/${study}" -name "$name_gwas_munge_file")
    match_file=${matches[0]}
    if [[ -n $matches && -f $match_file ]]; then
      secondaries+=("$study;${path_gwas_source_container}/${study}/${name_gwas_munge_file}")
    fi
  fi
done

# Assemble array of batch instance details.
comparisons=()
for primary in "${primaries[@]}"; do
  for secondary in "${secondaries[@]}"; do
    comparisons+=("${primary};${secondary}")
  done
done

################################################################################
# Append custom comparisons that do not follow the same pattern.

##########
# cohorts_models_linear_measurement
if true; then
  # Females to Males.
  pairs+=("female_vitamin_d_log;male_vitamin_d_log")
  pairs+=("female_albumin;male_albumin")
  pairs+=("female_steroid_globulin_log;male_steroid_globulin_log")
  pairs+=("female_oestradiol_log;male_oestradiol_log")
  pairs+=("female_testosterone_log;male_testosterone_log")
  # Premenopause Females to Postmenopause Females.
  pairs+=("female_premenopause_vitamin_d_log;female_postmenopause_vitamin_d_log")
  pairs+=("female_premenopause_albumin;female_postmenopause_albumin")
  pairs+=("female_premenopause_steroid_globulin_log;female_postmenopause_steroid_globulin_log")
  pairs+=("female_premenopause_oestradiol_log;female_postmenopause_oestradiol_log")
  pairs+=("female_premenopause_testosterone_log;female_postmenopause_testosterone_log")
  # Younger Males to Older Males.
  pairs+=("male_age_low_vitamin_d_log;male_age_high_vitamin_d_log")
  pairs+=("male_age_low_albumin;male_age_high_albumin")
  pairs+=("male_age_low_steroid_globulin_log;male_age_high_steroid_globulin_log")
  pairs+=("male_age_low_oestradiol_log;male_age_high_oestradiol_log")
  pairs+=("male_age_low_testosterone_log;male_age_high_testosterone_log")
fi

##########
# cohorts_models_linear_imputation
if false; then
  # Females to Males.
  pairs+=("female_vitamin_d_imputation_log;male_vitamin_d_imputation_log")
  pairs+=("female_albumin_imputation;male_albumin_imputation")
  pairs+=("female_steroid_globulin_imputation_log;male_steroid_globulin_imputation_log")
  pairs+=("female_oestradiol_imputation;male_oestradiol_imputation")
  pairs+=("female_testosterone_imputation;male_testosterone_imputation")
  # Premenopause Females to Postmenopause Females.
  pairs+=("female_premenopause_vitamin_d_imputation_log;female_postmenopause_vitamin_d_imputation_log")
  pairs+=("female_premenopause_albumin_imputation;female_postmenopause_albumin_imputation")
  pairs+=("female_premenopause_steroid_globulin_imputation_log;female_postmenopause_steroid_globulin_imputation_log")
  pairs+=("female_premenopause_oestradiol_imputation;female_postmenopause_oestradiol_imputation")
  pairs+=("female_premenopause_testosterone_imputation;female_postmenopause_testosterone_imputation")
  # Younger Males to Older Males.
  pairs+=("male_age_low_vitamin_d_imputation_log;male_age_high_vitamin_d_imputation_log")
  pairs+=("male_age_low_albumin_imputation;male_age_high_albumin_imputation")
  pairs+=("male_age_low_steroid_globulin_imputation_log;male_age_high_steroid_globulin_imputation_log")
  pairs+=("male_age_low_oestradiol_imputation;male_age_high_oestradiol_imputation")
  pairs+=("male_age_low_testosterone_imputation;male_age_high_testosterone_imputation")
fi

for pair in "${pairs[@]}"; do
  IFS=";" read -r -a array <<< "${pair}"
  study_primary="${array[0]}"
  study_secondary="${array[1]}"
  comparisons+=("${study_primary};${path_gwas_source_container}/${study_primary}/${name_gwas_munge_file};${study_secondary};${path_gwas_source_container}/${study_secondary}/${name_gwas_munge_file}")
done

################################################################################
# Drive genetic correlations across comparisons.
# Format for array of comparisons.
# "study_primary;path_gwas_primary_munge_suffix;study_secondary;path_gwas_secondary_munge_suffix"

for comparison in "${comparisons[@]}"; do

  ##############################################################################
  # Extract details for comparison.
  IFS=";" read -r -a array <<< "${comparison}"
  study_primary="${array[0]}"
  path_gwas_primary_munge_suffix="${array[1]}"
  study_secondary="${array[2]}"
  path_gwas_secondary_munge_suffix="${array[3]}"
  echo "----------"
  echo "primary study: ${study_primary}"
  echo "path: ${path_gwas_primary_munge_suffix}"
  echo "secondary study: ${study_secondary}"
  echo "path: ${path_gwas_secondary_munge_suffix}"

  if true; then
    ##############################################################################
    # LDSC Genetic Correlation.
    # Paths.
    path_genetic_reference="${path_dock}/access/genetic_reference"
    #study_primary=""
    #study_secondary=""
    #path_gwas_primary_munge_suffix=""
    #path_gwas_secondary_munge_suffix=""
    path_genetic_correlation_parent="${path_dock}/genetic_correlation/${cohorts_models}/${study_primary}/${study_secondary}"
    rm -r $path_genetic_correlation_parent
    mkdir -p $path_genetic_correlation_parent
    # Scripts.
    path_promiscuity_scripts="${path_process}/promiscuity/scripts"
    path_scripts_gwas_process="${path_promiscuity_scripts}/gwas_process"
    path_script_drive_ldsc_gwas_genetic_correlation="${path_scripts_gwas_process}/drive_ldsc_gwas_genetic_correlation.sh"
    # Parameters.
    report="true" # "true" or "false"
    /usr/bin/bash "$path_script_drive_ldsc_gwas_genetic_correlation" \
    $path_gwas_primary_munge_suffix \
    $path_gwas_secondary_munge_suffix \
    $path_genetic_correlation_parent \
    $path_genetic_reference \
    $report
  fi
done
