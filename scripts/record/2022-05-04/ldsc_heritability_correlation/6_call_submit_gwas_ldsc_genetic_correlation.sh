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

#cohorts_models="albumin_linear_2"                      #  22 Comparisons;  2 GWAS; 13 April 2022
#cohorts_models="albumin_linear_1"                      # 198 Comparisons; 18 GWAS; 13 April 2022
#cohorts_models="steroid_globulin_linear_2"             #  22 Comparisons;  2 GWAS; 13 April 2022
#cohorts_models="steroid_globulin_linear_1"             # 198 Comparisons; 18 GWAS; 13 April 2022
#cohorts_models="oestradiol_logistic"                   # __ Comparisons; __ GWAS; incomplete
#cohorts_models="oestradiol_linear_1"                   # 198 Comparisons; 18 GWAS; 13 April 2022
#cohorts_models="oestradiol_linear_2"                   # 132 Comparisons; 12 GWAS; 13 April 2022
#cohorts_models="oestradiol_bioavailable_linear"        # 198 Comparisons; 18 GWAS; 13 April 2022
#cohorts_models="oestradiol_free_linear"                # 198 Comparisons; 18 GWAS; 13 April 2022
#cohorts_models="testosterone_logistic"                 # __ Comparisons; __ GWAS; incomplete
#cohorts_models="testosterone_linear"                   # __ Comparisons; 54 GWAS; __ April 2022 <-- should be 594
#cohorts_models="testosterone_bioavailable_linear"      # 198 Comparisons; 18 GWAS; 13 April 2022
cohorts_models="testosterone_free_linear"              # __ Comparisons; 18 GWAS; __ April 2022

################################################################################
# Organize paths.
# Read private, local file paths.
cd ~/paths
path_process=$(<"./process_psychiatric_metabolism.txt")
path_dock="$path_process/dock"

path_genetic_reference="${path_dock}/access/genetic_reference"
path_genetic_correlation_container="${path_dock}/genetic_correlation"
name_gwas_munge_file="gwas_munge.sumstats.gz"
path_primary_gwas_munge_container="${path_dock}/gwas_ldsc_format_munge"
path_secondary_gwas_munge_container="${path_dock}/gwas_ldsc_munge/${cohorts_models}"
path_scripts_record="${path_process}/psychiatric_metabolism/scripts/record/2022-05-04/ldsc_heritability_correlation"

# Parameters.
report="true" # "true" or "false"

mkdir -p $path_genetic_correlation_container

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

###primaries+=("34255042_schmitz_2021_female;${path_primary_gwas_munge_container}/34255042_schmitz_2021_female/${name_gwas_munge_file}")
###primaries+=("34255042_schmitz_2021_male;${path_primary_gwas_munge_container}/34255042_schmitz_2021_male/${name_gwas_munge_file}")

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
if true; then
  # Assemble array of batch instance details.
  comparison_container="primary_versus_${cohorts_models}"
  for primary in "${primaries[@]}"; do
    for secondary in "${secondaries[@]}"; do
      comparisons+=("${comparison_container};${primary};${secondary}")
    done
  done
fi

##########
# Study pairs with special secondary studies for comparison to all primary studies.
if false; then
  # Collect special secondary studies.
  secondaries_special=()
  secondaries_special+=("34255042_schmitz_2021_female;${path_primary_gwas_munge_container}/34255042_schmitz_2021_female/${name_gwas_munge_file}")
  secondaries_special+=("34255042_schmitz_2021_male;${path_primary_gwas_munge_container}/34255042_schmitz_2021_male/${name_gwas_munge_file}")
  # Assemble array of batch instance details.
  comparison_container="primary_versus_34255042_schmitz_2021"
  for primary in "${primaries[@]}"; do
    for secondary_special in "${secondaries_special[@]}"; do
      comparisons+=("${comparison_container};${primary};${secondary_special}")
    done
  done
fi


##########
# Study pairs within the same container (path_primary_gwas_munge_container).
if false; then
  # Signal transformation.
  pairs+=("34255042_schmitz_2021_female;34255042_schmitz_2021_male")

  # Assemble array of batch instance details.
  comparison_container="34255042_schmitz_2021_female_versus_male"
  for pair in "${pairs[@]}"; do
    IFS=";" read -r -a array <<< "${pair}"
    study_primary="${array[0]}"
    study_secondary="${array[1]}"
    comparisons+=("${comparison_container};${study_primary};${path_primary_gwas_munge_container}/${study_primary}/${name_gwas_munge_file};${study_secondary};${path_primary_gwas_munge_container}/${study_secondary}/${name_gwas_munge_file}")
  done
fi

# TODO: need to update path to GWAS results for oestradiol_logistic in females and males...

##########
# Study pairs within different containers.
if false; then
  # Females to Males.
  pairs+=("34255042_schmitz_2021_female;${path_primary_gwas_munge_container}/34255042_schmitz_2021_female;male_vitamin_d;${path_dock}/gwas_ldsc_munge/oestradiol_logistic_long/female...")
  pairs+=("34255042_schmitz_2021_male;${path_primary_gwas_munge_container}/34255042_schmitz_2021_male;male_vitamin_d_log;${path_dock}/gwas_ldsc_munge/oestradiol_logistic_long/male...")

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

report="true"

if true; then
  # Execute comparisons in parallel compute batch job.
  # Initialize batch instances.
  path_batch_instances="${path_genetic_correlation_container}/batch_instances_${cohorts_models}.txt"
  rm $path_batch_instances
  # Write out batch instances.
  for comparison in "${comparisons[@]}"; do
    echo $comparison >> $path_batch_instances
  done
  # Read batch instances.
  readarray -t batch_instances < $path_batch_instances
  batch_instances_count=${#batch_instances[@]}
  echo "----------"
  echo "count of batch instances: " $batch_instances_count
  echo "first batch instance: " ${batch_instances[0]} # notice base-zero indexing
  echo "last batch instance: " ${batch_instances[$batch_instances_count - 1]}
  # Submit array batch to Sun Grid Engine.
  # Array batch indices must start at one (not zero).
  qsub -t 1-${batch_instances_count}:1 -o \
  "${path_genetic_correlation_container}/batch_${cohorts_models}_out.txt" -e "${path_genetic_correlation_container}/batch_${cohorts_models}_error.txt" \
  "${path_scripts_record}/7_run_batch_jobs_gwas_ldsc_genetic_correlation.sh" \
  $path_batch_instances \
  $batch_instances_count \
  $path_genetic_correlation_container \
  $path_genetic_reference \
  $path_scripts_record \
  $path_process \
  $report
else
  # Execute comparisons by simple, sequential iteration.
  for comparison_instance in "${comparisons[@]}"; do
    # Call driver script.
    /usr/bin/bash "${path_scripts_record}/8_run_drive_gwas_ldsc_genetic_correlation.sh" \
    $comparison_instance \
    $path_genetic_correlation_container \
    $path_genetic_reference \
    $path_process \
    $report
  done
fi



#
