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

cohorts_models="vitamin_d_linear"                   # 4 GWAS; GWAS job 3155905, status: complete; TCW started at ___ on 03 March 2022;

###cohorts_models="oestradiol_logistic"              # 24 GWAS; GWAS job 3202343, status: running; <-- priority!!!
#cohorts_models="oestradiol_bioavailable_linear"     # 18 GWAS; GWAS job 3149651, status: complete; TCW started at ___ on 03 March 2022;
#cohorts_models="oestradiol_free_linear"             # 18 GWAS; GWAS job 3149652, status: complete; TCW started at ___ on 03 March 2022;

###cohorts_models="testosterone_logistic"            # 24 GWAS; GWAS job 3202423, status: in queue;
#cohorts_models="testosterone_linear"                # 24 GWAS; GWAS job 3109689, status: complete; TCW started at ___ on 03 March 2022;
#cohorts_models="testosterone_bioavailable_linear"   # 18 GWAS; GWAS job 3149548, status: complete; TCW started at ___ on 03 March 2022;
#cohorts_models="testosterone_free_linear"           # 18 GWAS; GWAS job 3149549, status: complete; TCW started at ___ on 03 March 2022;

#cohorts_models="steroid_globulin_linear"            # 4 GWAS;  GWAS job 3155785, status: complete; TCW started at ___ on 03 March 2022;
###cohorts_models="steroid_globulin_sex_linear"      # 24 GWAS; GWAS job 3202509, status: in queue;
#cohorts_models="albumin_linear"                     # 4 GWAS;  GWAS job 3155786, status: complete; TCW started at ___ on 03 March 2022;

################################################################################
# Organize paths.
# Read private, local file paths.
cd ~/paths
path_process=$(<"./process_psychiatric_metabolism.txt")
path_dock="$path_process/dock"

path_genetic_correlation_container="${path_dock}/genetic_correlation"
name_gwas_munge_file="gwas_munge.sumstats.gz"
path_primary_gwas_munge_container="${path_dock}/gwas_ldsc_format_munge"
path_secondary_gwas_munge_container="${path_dock}/gwas_ldsc_munge/${cohorts_models}"

###########################################################################
# Define main comparisons.

# Define array of primary studies.
primaries=()
primaries+=("30482948_walters_2018_eur_unrel_meta;${path_primary_gwas_munge_container}/30482948_walters_2018_eur_unrel_meta/${name_gwas_munge_file}")
primaries+=("30482948_walters_2018_eur_unrel_genotype;${path_primary_gwas_munge_container}/30482948_walters_2018_eur_unrel_genotype/${name_gwas_munge_file}")
primaries+=("30482948_walters_2018_female;${path_primary_gwas_munge_container}/30482948_walters_2018_female/${name_gwas_munge_file}")
primaries+=("30482948_walters_2018_male;${path_primary_gwas_munge_container}/30482948_walters_2018_male/${name_gwas_munge_file}")
primaries+=("30643251_liu_2019_alcohol_all;${path_primary_gwas_munge_container}/30643251_liu_2019_alcohol_all/${name_gwas_munge_file}")
primaries+=("30643251_liu_2019_alcohol_no_ukb;${path_primary_gwas_munge_container}/30643251_liu_2019_alcohol_no_ukb/${name_gwas_munge_file}")
primaries+=("30718901_howard_2019;${path_primary_gwas_munge_container}/30718901_howard_2019/${name_gwas_munge_file}")
primaries+=("34002096_mullins_2021_all;${path_primary_gwas_munge_container}/34002096_mullins_2021_all/${name_gwas_munge_file}")
primaries+=("34002096_mullins_2021_bpd1;${path_primary_gwas_munge_container}/34002096_mullins_2021_bpd1/${name_gwas_munge_file}")
primaries+=("34002096_mullins_2021_bpd2;${path_primary_gwas_munge_container}/34002096_mullins_2021_bpd2/${name_gwas_munge_file}")
primaries+=("00000000_ripke_2022;${path_primary_gwas_munge_container}/00000000_ripke_2022/${name_gwas_munge_file}")

#primaries+=("34255042_schmitz_2021_female;${path_primary_gwas_munge_container}/34255042_schmitz_2021_female/${name_gwas_munge_file}")
#primaries+=("34255042_schmitz_2021_male;${path_primary_gwas_munge_container}/34255042_schmitz_2021_male/${name_gwas_munge_file}")

# Define array of secondary studies.
secondaries=()
# Iterate on directories for GWAS on cohorts and hormones.
cd $path_secondary_gwas_munge_container
for path_directory in `find . -maxdepth 1 -mindepth 1 -type d -not -name .`; do
  if [ -d "$path_directory" ]; then
    # Current content item is a directory.
    # Extract directory's base name.
    study="$(basename -- $path_directory)"
    #echo $directory
    # Determine whether directory contains valid GWAS summary statistics.
    matches=$(find "${path_secondary_gwas_munge_container}/${study}" -name "$name_gwas_munge_file")
    match_file=${matches[0]}
    if [[ -n $matches && -f $match_file ]]; then
      secondaries+=("$study;${path_secondary_gwas_munge_container}/${study}/${name_gwas_munge_file}")
    fi
  fi
done

################################################################################
# Collect comparisons between studies.

# Initialize array of comparisons.
comparisons=()

##########
# Comparison pairs of secondary studies for comparison to all primary studies.
if false; then
  # Assemble array of batch instance details.
  comparison_container="${cohorts_models}_primary_secondary"
  for primary in "${primaries[@]}"; do
    for secondary in "${secondaries[@]}"; do
      comparisons+=("${comparison_container};${primary};${secondary}")
    done
  done
fi

##########
# Study pairs with special secondary studies for comparison to all primary studies.
if true; then
  # Collect special secondary studies.
  secondaries_special=()
  secondaries_special+=("34255042_schmitz_2021_female;${path_primary_gwas_munge_container}/34255042_schmitz_2021_female/${name_gwas_munge_file}")
  secondaries_special+=("34255042_schmitz_2021_male;${path_primary_gwas_munge_container}/34255042_schmitz_2021_male/${name_gwas_munge_file}")
  # Assemble array of batch instance details.
  comparison_container="comparisons_primary_34255042_schmitz_2021"
  for primary in "${primaries[@]}"; do
    for secondary_special in "${secondaries_special[@]}"; do
      comparisons+=("${comparison_container};${primary};${secondary_special}")
    done
  done
fi


##########
# Study pairs within the same container (white_female_male_priority_male_linear).
if false; then
  # Signal transformation.
  pairs+=("joint_vitamin_d;joint_vitamin_d_log")
  pairs+=("joint_vitamin_d_imputation;joint_vitamin_d_imputation_log")
  # Model.
  pairs+=("unadjust_vitamin_d_imputation_log;joint_vitamin_d_imputation_log")

  # Assemble array of batch instance details.
  comparison_container="${cohorts_models}_secondary_pairs_model"
  for pair in "${pairs[@]}"; do
    IFS=";" read -r -a array <<< "${pair}"
    study_primary="${array[0]}"
    study_secondary="${array[1]}"
    comparisons+=("${comparison_container};${study_primary};${path_secondary_gwas_munge_container}/${study_primary}/${name_gwas_munge_file};${study_secondary};${path_secondary_gwas_munge_container}/${study_secondary}/${name_gwas_munge_file}")
  done
fi

##########
# Study pairs within different containers.
if false; then
  # Females to Males.
  pairs+=("female_vitamin_d;${path_dock}/gwas_ldsc_munge/white_female_linear/joint_vitamin_d;male_vitamin_d;${path_dock}/gwas_ldsc_munge/white_male_linear/joint_vitamin_d")
  pairs+=("female_vitamin_d_log;${path_dock}/gwas_ldsc_munge/white_female_linear/joint_vitamin_d_log;male_vitamin_d_log;${path_dock}/gwas_ldsc_munge/white_male_linear/joint_vitamin_d_log")
  pairs+=("female_vitamin_d_imputation;${path_dock}/gwas_ldsc_munge/white_female_linear/joint_vitamin_d_imputation;male_vitamin_d_imputation;${path_dock}/gwas_ldsc_munge/white_male_linear/joint_vitamin_d_imputation")
  pairs+=("female_vitamin_d_imputation_log;${path_dock}/gwas_ldsc_munge/white_female_linear/joint_vitamin_d_imputation_log;male_vitamin_d_imputation_log;${path_dock}/gwas_ldsc_munge/white_male_linear/joint_vitamin_d_imputation_log")

  # Assemble array of batch instance details.
  comparison_container="white_secondary_pairs_female_male"
  for pair in "${pairs[@]}"; do
    IFS=";" read -r -a array <<< "${pair}"
    study_primary="${array[0]}"
    path_primary="${array[1]}"
    study_secondary="${array[2]}"
    path_secondary="${array[3]}"
    comparisons+=("${comparison_container};${study_primary};${path_primary}/${name_gwas_munge_file};${study_secondary};${path_secondary}/${name_gwas_munge_file}")
  done
fi


################################################################################
# Drive genetic correlations across comparisons.
# Format for array of comparisons.
# "study_primary;path_gwas_primary_munge_suffix;study_secondary;path_gwas_secondary_munge_suffix"

for comparison in "${comparisons[@]}"; do

  ##############################################################################
  # Extract details for comparison.
  IFS=";" read -r -a array <<< "${comparison}"
  comparison_container="${array[0]}"
  study_primary="${array[1]}"
  path_gwas_primary_munge_suffix="${array[2]}"
  study_secondary="${array[3]}"
  path_gwas_secondary_munge_suffix="${array[4]}"
  echo "----------"
  echo "comparison container: ${comparison_container}"
  echo "primary study: ${study_primary}"
  echo "path: ${path_gwas_primary_munge_suffix}"
  echo "secondary study: ${study_secondary}"
  echo "path: ${path_gwas_secondary_munge_suffix}"

  if true; then
    ##############################################################################
    # LDSC Genetic Correlation.
    # Paths.
    path_genetic_reference="${path_dock}/access/genetic_reference"
    path_genetic_correlation_parent="${path_genetic_correlation_container}/${comparison_container}/${study_primary}/${study_secondary}"
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
