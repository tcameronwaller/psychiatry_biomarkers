#!/bin/bash

###########################################################################
###########################################################################
###########################################################################
# ...
###########################################################################
###########################################################################
###########################################################################

################################################################################
# Activate Virtual Environment.

# Read private, local file paths.
#echo "read private file path variables and organize paths..."
cd ~/paths
path_tools=$(<"./waller_tools.txt")
path_environment_ldsc="${path_tools}/python/environments/ldsc"
source "${path_environment_ldsc}/bin/activate"
which python2
sleep 5s

################################################################################
# Organize paths.
# Read private, local file paths.
cd ~/paths
path_ldsc=$(<"./tools_ldsc.txt")
path_process=$(<"./process_psychiatric_metabolism.txt")

path_dock="$path_process/dock"
path_genetic_reference="${path_dock}/access/genetic_reference"
path_gwas="${path_dock}/gwas"
path_gwas_parent="${path_gwas}/cohorts_models" # selection
path_heritability="${path_dock}/heritability"
path_heritability_parent="${path_heritability}/cohorts_models" # selection

path_scripts_record="$path_process/psychiatric_metabolism/scripts/record/2021-07-02"
path_promiscuity_scripts="${path_process}/promiscuity/scripts"

# Initialize directories.
rm -r $path_heritability_parent
mkdir -p $path_heritability_parent

###########################################################################
# Define explicit inclusions and exclusions.
# Use inclusions to run procedure for a few specific cohort-hormone combinations that are missing from the set.
# Use exclusions to omit a few cohort-hormone combinations that are not complete yet.

delimiter=" "
IFS=${delimiter}
exclusions=()
#exclusions+=("female_combination_unadjust_albumin_log")
unset IFS

###########################################################################
# Execute procedure.

# Initialize batch instances.
path_batch_instances="${path_gwas_parent}/batch_instances.txt"
rm $path_batch_instances

# Iterate on directories for GWAS on cohorts and hormones.
pattern_gwas_check_file="report.*.glm.linear" # do not expand with full path yet
cd $path_gwas_parent
for path_directory in `find . -maxdepth 1 -mindepth 1 -type d -not -name .`; do
  if [ -d "$path_directory" ]; then
    # Current content item is a directory.
    directory="$(basename -- $path_directory)"

    # Determine specific inclusions or exclusions.
    # inclusions: [[ " ${inclusions[@]} " =~ "${directory}" ]]
    # exclusions: [[ ! " ${exclusions[@]} " =~ "${directory}" ]]
    if [[ ! " ${exclusions[@]} " =~ "${directory}" ]]; then
      # Determine whether directory contains valid GWAS summary statistics
      # across chromosomes.
      # Check for chromosome 22, assuming that all chromosomes completed
      # sequentially.
      matches=("${path_gwas_parent}/${directory}/chromosome_22/${pattern_gwas_check_file}")
      path_gwas_check_file="${matches[0]}"
      echo $path_gwas_check_file
      echo $directory
      echo $directory >> $path_batch_instances
      #if [[ -f "$path_gwas_check_file" ]]; then
      #  echo $directory
      #  echo $directory >> $path_batch_instances
      #fi
    fi
  fi
done

# Read batch instances.
readarray -t batch_instances < $path_batch_instances
batch_instances_count=${#batch_instances[@]}
echo "----------"
echo "count of batch instances: " $batch_instances_count
echo "first batch instance: " ${batch_instances[0]} # notice base-zero indexing
echo "last batch instance: " ${batch_instances[batch_instances_count - 1]}

# Execute batch with grid scheduler.
if false; then
  report="true" # "true" or "false"
  # Submit array batch to Sun Grid Engine.
  # Array batch indices must start at one (not zero).
  qsub -t 1-${batch_instances_count}:1 -o \
  "${path_gwas_parent}/out.txt" -e "${path_gwas_parent}/error.txt" \
  "${path_promiscuity_scripts}/run_plink_gwas_collection_format_munge_heritability.sh" \
  $path_batch_instances \
  $batch_instances_count \
  $path_gwas_parent \
  $path_heritability_parent \
  $path_genetic_reference \
  $path_promiscuity_scripts \
  $path_ldsc \
  $report
fi

# Execute iteratively on head node.
if true; then
  for directory in "${batch_instances[@]}"; do
    ###########################################################################
    # Organize variables.
    path_script_gwas_collect_concatenate="${path_promiscuity_scripts}/collect_concatenate_gwas_chromosomes.sh"
    path_promiscuity_scripts_ldsc_heritability="${path_promiscuity_scripts}/ldsc_genetic_heritability_correlation"
    path_script_format_munge_heritability="${path_promiscuity_scripts_ldsc_heritability}/format_munge_gwas_heritability_ldsc.sh"
    path_scripts_format="${path_promiscuity_scripts}/format_gwas_ldsc"
    path_script_gwas_format="${path_scripts_format}/format_gwas_ldsc_plink_linear.sh"

    ###########################################################################
    # Execute procedure.

    # Concatenate GWAS across chromosomes.
    if false; then
      # Organize variables.
      pattern_source_file="report.*.glm.linear" # do not expand with full path yet
      path_source_directory="${path_gwas_parent}/${directory}"
      chromosome_start=1
      chromosome_end=22
      path_gwas_concatenation="${path_source_directory}/gwas_concatenation.txt"
      path_gwas_concatenation_compress="${path_source_directory}/gwas_concatenation.txt.gz"
      /usr/bin/bash "$path_script_gwas_collect_concatenate" \
      $pattern_source_file \
      $path_source_directory \
      $chromosome_start \
      $chromosome_end \
      $path_gwas_concatenation \
      $path_gwas_concatenation_compress \
      $report
    fi

    # Format and munge GWAS summary statistics.
    # Estimate genotype heritability.
    if true; then
      # Organize paths.
      path_gwas_study="${path_gwas_parent}/${directory}"
      path_heritability_study="${path_heritability_parent}/${directory}"
      # Initialize directories.
      #mkdir -p $path_gwas_study
      mkdir -p $path_heritability_study
      # Organize variables.
      study="${directory}"
      name_prefix="null"
      path_source_file="${path_gwas_study}/gwas_concatenation.txt.gz"
      report="true" # "true" or "false"
      /usr/bin/bash "$path_script_format_munge_heritability" \
      $study \
      $name_prefix \
      $path_source_file \
      $path_genetic_reference \
      $path_gwas_study \
      $path_heritability_study \
      $path_script_gwas_format \
      $path_promiscuity_scripts \
      $path_ldsc \
      $report
    fi
  done
fi



################################################################################
# Deactivate Virtual Environment.

deactivate
which python2
