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
path_gwas_summaries=$(<"./gwas_summaries_waller_metabolism.txt")
path_process=$(<"./process_psychiatric_metabolism.txt")

path_dock="$path_process/dock"
path_genetic_reference="${path_dock}/access/genetic_reference"
path_gwas="${path_dock}/gwas"
path_heritability="${path_dock}/heritability"
#path_genetic_correlation="${path_dock}/genetic_correlation"

path_promiscuity_scripts="${path_process}/promiscuity/scripts"
path_promiscuity_scripts_ldsc_heritability="${path_promiscuity_scripts}/ldsc_genetic_heritability_correlation"
path_scripts_format="${path_promiscuity_scripts}/format_gwas_ldsc"
#path_scripts_record="$path_process/psychiatric_metabolism/scripts/record/2021-05-04"


################################################################################
# Organize variables.

metabolite_study="27005778_kettunen_2016"
path_metabolite_gwas_source_directory="${path_gwas_summaries}/${metabolite_study}"
metabolite_file_pattern="Summary_statistics_MAGNETIC_*.txt.gz" # do not expand with full path yet
metabolite_file_prefix="Summary_statistics_MAGNETIC_" # file name prefix before metabolite identifier or "null"
metabolite_file_suffix=".txt.gz" # file name suffix after metabolite identifier or "null"

path_study_gwas="${path_gwas}/${metabolite_study}"
path_batch_instances="${path_study_gwas}/batch_instances.txt"
path_study_heritability="${path_heritability}/${metabolite_study}"

path_script_gwas_format="${path_scripts_format}/format_gwas_ldsc_${metabolite_study}.sh"

# Initialize directories.
rm -r $path_study_gwas
rm -r $path_study_heritability
mkdir -p $path_study_gwas
mkdir -p $path_study_heritability

# Organize instances for iteration.
echo "----------------------------------------------------------------------"
echo "Organize array of batch instances."
echo "----------------------------------------------------------------------"
# Collect batch instances.
#cd $path_source
#metabolite_files=(metabolite_*_meta_analysis_gwas.csv.gz)
#for path_file in "${metabolite_files[@]}"; do
#    echo $path_file >> $path_metabolites/metabolite_files.txt
#done
# Define glob pattern for file paths.
# This definition expands to an array of all files in the path directory that
# match the pattern.
path_pattern="${path_metabolite_gwas_source_directory}/${metabolite_file_pattern}"
# Iterate on all files and directories in parent directory.
rm $path_batch_instances
for path_file in $path_metabolite_gwas_source_directory/*; do
  if [ -f "$path_file" ]; then
    # Current content item is a file.
    # Compare to glob pattern to recognize relevant files.
    #echo $path_file
    #echo $path_file >> $path_batch_instances
    if [[ "$path_file" == ${path_pattern} ]]; then
      # File name matches glob pattern.
      # Include full path to file in batch instances.
      #echo $path_file
      echo $path_file >> $path_batch_instances
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

if true; then
  for path_source_file in "${batch_instances[@]}"; do
    study=${metabolite_study}
    report="true" # "true" or "false"
    /usr/bin/bash "$path_promiscuity_scripts_ldsc_heritability/format_munge_gwas_heritability_ldsc.sh" \
    $study \
    $path_source_file \
    $path_genetic_reference \
    $path_gwas \
    $path_heritability \
    $path_script_gwas_format \
    $path_promiscuity_scripts \
    $path_ldsc \
    $report
  done
fi
if false; then
  # Submit array batch to Sun Grid Engine.
  # Array batch indices must start at one (not zero).
  echo "----------------------------------------------------------------------"
  echo "Submit array of batches to Sun Grid Engine."
  echo "----------------------------------------------------------------------"
  qsub -t 1-${batch_instances_count}:1 -o \
  "${path_study_gwas}/out.txt" -e "${path_study_gwas}/error.txt" \
  "${path_scripts_record}/6_run_batch_organize_gwas_ldsc_heritability.sh" \
  $path_batch_instances \
  $batch_instances_count \
  $phenotype_study \
  $metabolite_study \
  $name_prefix \
  $name_suffix \
  $path_genetic_reference \
  $path_phenotype_gwas \
  $path_study_gwas \
  $path_study_heritability \
  $path_study_genetic_correlation \
  $path_scripts_record \
  $path_script_gwas_format \
  $path_promiscuity_scripts \
  $report
fi
