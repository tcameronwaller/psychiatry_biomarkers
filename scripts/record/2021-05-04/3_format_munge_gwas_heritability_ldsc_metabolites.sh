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

path_promiscuity_scripts="${path_process}/promiscuity/scripts"
path_promiscuity_scripts_ldsc_heritability="${path_promiscuity_scripts}/ldsc_genetic_heritability_correlation"
path_format_munge_gwas_heritability_ldsc="${path_promiscuity_scripts_ldsc_heritability}/format_munge_gwas_heritability_ldsc.sh"
path_scripts_format="${path_promiscuity_scripts}/format_gwas_ldsc"
path_scripts_record="$path_process/psychiatric_metabolism/scripts/record/2021-05-04"


################################################################################
# Organize variables.

if false; then
  metabolite_study="24816252_shin_2014"
  path_metabolite_gwas_source_directory="${path_gwas_summaries}/${metabolite_study}/metabolites_meta"
  metabolite_file_pattern="*.metal.pos.txt.gz" # do not expand with full path yet
  metabolite_file_prefix="null" # file name prefix before metabolite identifier or "null"
  metabolite_file_suffix=".metal.pos.txt.gz" # file name suffix after metabolite identifier or "null"
fi

if true; then
  metabolite_study="27005778_kettunen_2016"
  path_metabolite_gwas_source_directory="${path_gwas_summaries}/${metabolite_study}"
  metabolite_file_pattern="Summary_statistics_MAGNETIC_*.txt.gz" # do not expand with full path yet
  metabolite_file_prefix="Summary_statistics_MAGNETIC_" # file name prefix before metabolite identifier or "null"
  metabolite_file_suffix=".txt.gz" # file name suffix after metabolite identifier or "null"
fi

if false; then
  metabolite_study="33437055_panyard_2021"
  path_metabolite_gwas_source_directory="${path_gwas_summaries}/${metabolite_study}"
  metabolite_file_pattern="metabolite_*_meta_analysis_gwas.csv.gz" # do not expand with full path yet
  metabolite_file_prefix="metabolite_" # file name prefix before metabolite identifier or "null"
  metabolite_file_suffix="_meta_analysis_gwas.csv.gz" # file name suffix after metabolite identifier or "null"
fi

path_gwas_study="${path_gwas}/${metabolite_study}"
path_batch_instances="${path_gwas_study}/batch_instances.txt"
path_heritability_study="${path_heritability}/${metabolite_study}"

path_script_gwas_format="${path_scripts_format}/format_gwas_ldsc_${metabolite_study}.sh"

# Initialize directories.
rm -r $path_gwas_study
rm -r $path_heritability_study
mkdir -p $path_gwas_study
mkdir -p $path_heritability_study

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

report="true" # "true" or "false"

# Execute batch with grid scheduler.
if false; then
  # Submit array batch to Sun Grid Engine.
  # Array batch indices must start at one (not zero).
  qsub -t 1-${batch_instances_count}:1 -o \
  "${path_gwas_study}/out.txt" -e "${path_gwas_study}/error.txt" \
  "${path_scripts_record}/4_run_batch_metabolite_format_munge_gwas_heritability_ldsc.sh" \
  $path_batch_instances \
  $batch_instances_count \
  $metabolite_study \
  $metabolite_file_prefix \
  $metabolite_file_suffix \
  $path_genetic_reference \
  $path_gwas_study \
  $path_heritability_study \
  $path_script_gwas_format \
  $path_format_munge_gwas_heritability_ldsc \
  $path_promiscuity_scripts \
  $path_ldsc \
  $report
fi

# Execute batch iteratively without grid scheduler.
if false; then
  for path_source_file in "${batch_instances[@]}"; do
    # Determine file name.
    #file_name=$path_source_file
    file_name="$(basename -- $path_source_file)"
    # Determine metabolite identifier.
    # Refer to documnetation for test: https://www.freebsd.org/cgi/man.cgi?test
    # Bash script more or less ignores empty string argument, so it does not
    # work well to pass an empty string as an argument.
    # Instead use a non-empty string, such as "null".
    # if [[ ! -z "$name_prefix" ]]; then
    metabolite=${file_name}
    if [[ "$metabolite_file_prefix" != "null" ]]; then
      metabolite=${metabolite/$metabolite_file_prefix/""}
    fi
    if [[ "$metabolite_file_suffix" != "null" ]]; then
      metabolite=${metabolite/$metabolite_file_suffix/""}
    fi
    # Report.
    if [[ "$report" == "true" ]]; then
      echo "----------------------------------------------------------------------"
      echo "----------------------------------------------------------------------"
      echo "----------------------------------------------------------------------"
      echo "metabolite study: " $metabolite_study
      echo "path to metabolite file: " $path_source_file
      echo "file: " $file_name
      echo "metabolite: " $metabolite
      echo "----------"
    fi

    study=${metabolite_study}
    name_prefix=${metabolite} # file name prefix or "null"
    /usr/bin/bash "$path_format_munge_gwas_heritability_ldsc" \
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
  done
fi
