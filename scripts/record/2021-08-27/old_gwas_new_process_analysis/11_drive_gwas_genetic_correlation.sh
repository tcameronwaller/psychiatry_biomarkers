#!/bin/bash

###########################################################################
###########################################################################
###########################################################################
# ...
###########################################################################
###########################################################################
###########################################################################

################################################################################
# Organize paths.
# Read private, local file paths.
cd ~/paths
path_process=$(<"./process_psychiatric_metabolism.txt")
path_dock="$path_process/dock"

###########################################################################
# Define main comparisons.

# Define array of primary studies.
primaries=()
primaries+=("30124842_yengo_2018;${path_dock}/gwas_process/30124842_yengo_2018/gwas_munge.sumstats.gz")

# Define array of secondary studies.
secondaries=()
# Iterate on directories for GWAS on cohorts and hormones.
path_gwas_source_container="${path_dock}/gwas_process/cohorts_models"
name_gwas_munge_file="gwas_munge.sumstats.gz"
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
# Define custom comparisons.

# TODO: see script 8 from 11 June 2021

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
    path_genetic_correlation_parent="${path_dock}/genetic_correlation/${study_primary}/cohorts_models/${study_secondary}"
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
