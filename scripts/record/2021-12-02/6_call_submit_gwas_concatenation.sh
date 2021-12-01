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

cohorts_models="body_white_bipolar_strict"          # 12 GWAS; TCW started at ___ on 30 November 2021
#cohorts_models="body_white_bipolar_loose"          # 12 GWAS; TCW started at ___ on 30 November 2021

pattern_gwas_report_file="report.*.glm.linear" # do not expand with full path yet
chromosome_x="false" # whether to collect GWAS summary statistics report for Chromosome X

################################################################################
# Organize paths.
# Read private, local file paths.
cd ~/paths
path_process=$(<"./process_psychiatric_metabolism.txt")
path_dock="$path_process/dock"

path_gwas_source_container="${path_dock}/gwas_raw/${cohorts_models}"
path_gwas_target_container="${path_dock}/gwas_concatenation_no_x/${cohorts_models}"

path_scripts_record="$path_process/psychiatric_metabolism/scripts/record/2021-12-02"
path_batch_instances="${path_gwas_target_container}/batch_instances_concatenation.txt"

###########################################################################
# Execute procedure.

# Initialize batch instances.
rm -r $path_gwas_target_container
mkdir -p $path_gwas_target_container
rm $path_batch_instances

# Iterate on directories for GWAS on cohorts and hormones.
pattern_gwas_concatenation_file="gwas_concatenation.txt.gz"
cd $path_gwas_source_container
for path_directory in `find . -maxdepth 1 -mindepth 1 -type d -not -name .`; do
  if [ -d "$path_directory" ]; then
    # Current content item is a directory.
    # Extract directory's base name.
    study="$(basename -- $path_directory)"
    #echo $directory
    # Determine whether directory contains valid GWAS summary statistics
    # across chromosomes.
    # Check for chromosome 22, assuming that all chromosomes completed
    # sequentially.
    matches_chromosome=$(find "${path_gwas_source_container}/${study}/chromosome_22" -name "$pattern_gwas_report_file")
    match_chromosome_file=${matches_chromosome[0]}
    if [[ -n $matches_chromosome && -f $match_chromosome_file ]]; then
      echo "----------"
      echo "----------"
      echo "----------"
      echo "Found chromosome-22 sum stats for: ${study}"
      echo $matches_chromosome
      echo "Found match file: ${match_chromosome_file}"
      # Determine whether target directory already contains a concatenation file.
      mkdir -p "${path_gwas_target_container}/${study}"
      matches_concatenation=$(find "${path_gwas_target_container}/${study}" -name "$pattern_gwas_concatenation_file")
      match_concatenation_file=${matches_concatenation[0]}
      if [[ -n $matches_concatenation && -f $match_concatenation_file ]]; then
        echo "-----"
        echo "Found GWAS concatenation for: ${study}"
        echo "Found match file: ${match_concatenation_file}"
      else
        echo "Study has chromosome GWAS summary statistics but no concatenation."
        echo $study >> $path_batch_instances
      fi
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
  "${path_gwas_target_container}/out_concatenation.txt" -e "${path_gwas_target_container}/error_concatenation.txt" \
  "${path_scripts_record}/7_run_batch_jobs_gwas_concatenation.sh" \
  $path_batch_instances \
  $batch_instances_count \
  $pattern_gwas_report_file \
  $path_gwas_source_container \
  $path_gwas_target_container \
  $path_scripts_record \
  $path_process \
  $chromosome_x
fi
