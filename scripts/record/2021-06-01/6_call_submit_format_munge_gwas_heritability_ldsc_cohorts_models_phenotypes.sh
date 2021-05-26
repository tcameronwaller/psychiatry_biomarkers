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
path_ldsc=$(<"./tools_ldsc.txt")
path_process=$(<"./process_psychiatric_metabolism.txt")

path_dock="$path_process/dock"
path_genetic_reference="${path_dock}/access/genetic_reference"
path_gwas="${path_dock}/gwas"
path_gwas_cohorts_models="${path_gwas}/cohorts_models"
path_heritability="${path_dock}/heritability"
path_heritability_cohorts_models="${path_heritability}/cohorts_models"

path_scripts_record="$path_process/psychiatric_metabolism/scripts/record/2021-06-01"
path_promiscuity_scripts="${path_process}/promiscuity/scripts"
path_script_gwas_collect_concatenate="${path_promiscuity_scripts}/collect_concatenate_gwas_chromosomes.sh"
path_promiscuity_scripts_ldsc_heritability="${path_promiscuity_scripts}/ldsc_genetic_heritability_correlation"
path_script_format_munge_heritability="${path_promiscuity_scripts_ldsc_heritability}/format_munge_gwas_heritability_ldsc.sh"
path_scripts_format="${path_promiscuity_scripts}/format_gwas_ldsc"
path_script_gwas_format="${path_scripts_format}/format_gwas_ldsc_plink_linear.sh"

# Initialize directories.
rm -r $path_heritability_cohorts_models
mkdir -p $path_heritability_cohorts_models

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
path_batch_instances="${path_gwas_cohorts_models}/batch_instances.txt"
rm $path_batch_instances

# Iterate on directories for GWAS on cohorts and hormones.
cd $path_gwas_cohorts_models
for path_directory in `find . -maxdepth 1 -mindepth 1 -type d -not -name .`; do
  if [ -d "$path_directory" ]; then
    # Current content item is a directory.
    directory="$(basename -- $path_directory)"

    # Determine specific inclusions or exclusions.
    # inclusions: [[ " ${inclusions[@]} " =~ "${directory}" ]]
    # exclusions: [[ ! " ${exclusions[@]} " =~ "${directory}" ]]
    if [[ ! " ${exclusions[@]} " =~ "${directory}" ]]; then

      echo $directory
      echo $directory >> $path_batch_instances

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
if true; then

  report="true" # "true" or "false"
  # Submit array batch to Sun Grid Engine.
  # Array batch indices must start at one (not zero).
  qsub -t 1-${batch_instances_count}:1 -o \
  "${path_gwas_cohorts_models}/out.txt" -e "${path_gwas_cohorts_models}/error.txt" \
  "${path_scripts_record}/7_format_munge_gwas_heritability_ldsc_cohorts_models_phenotypes.sh" \
  $path_batch_instances \
  $batch_instances_count \
  $path_gwas_cohorts_models \
  $path_heritability_cohorts_models \
  $path_genetic_reference \
  $path_script_gwas_collect_concatenate \
  $path_script_gwas_format \
  $path_script_format_munge_heritability \
  $path_promiscuity_scripts \
  $path_ldsc \
  $report
fi
