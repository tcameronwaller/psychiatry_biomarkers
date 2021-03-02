#!/bin/bash

###########################################################################
###########################################################################
###########################################################################
# ...
###########################################################################
###########################################################################
###########################################################################

# Organize variables.
path_source=$1 # full path to source directory of GWAS summary statistics
file_pattern=$2 # glob pattern by which to recognize relevant files in source directory
path_destination_parent=$3 # full path to destination directory
path_script_gwas_organization=$4 # full path to script to use for format organization
path_scripts=$5 # full path to scripts for current implementation pipeline
path_promiscuity_scripts=$6 # full path to scripts from promiscuity package

# Initialize directories.
rm -r $path_destination_parent
if [ ! -d $path_destination_parent ]; then
    # Directory does not already exist.
    # Create directory.
    mkdir -p $path_destination_parent
fi

# Organize instances for iteration.
echo "----------------------------------------------------------------------"
echo "Organize array of batch instances."
echo "----------------------------------------------------------------------"
# Collect batch instances.
path_batch_instances="${path_destination_parent}/batch_instances.txt"
#cd $path_source
#metabolite_files=(metabolite_*_meta_analysis_gwas.csv.gz)
#for path_file in "${metabolite_files[@]}"; do
#    echo $path_file >> $path_metabolites/metabolite_files.txt
#done
# Iterate on all files and directories in parent directory.
for path_file in $path_source/*; do
  if [ -f "$path_file" ]; then
    # Current content item is a file.
    # Compare to glob pattern to recognize relevant files.
    echo $path_file >> $path_batch_instances
    if [[ "$path_file" == ${file_pattern} ]]; then
      # File name matches glob pattern.
      # Include full path to file in batch instances.
      echo $path_file >> $path_batch_instances
    fi
  fi
done

# Read batch instances.
readarray -t batch_instances < $path_batch_instances
batch_instances_count=${#batch_instances[@]}
echo "----------"
echo "count of batch instances: " $batch_instances_count
# Submit array batch to Sun Grid Engine.
# Array batch indices cannot start at zero.
echo "----------------------------------------------------------------------"
echo "Submit array of batches to Sun Grid Engine."
echo "----------------------------------------------------------------------"
if false; then
  qsub -t 1-${batch_instances_count}:1 -o \
  "$path_destination_parent/out.txt" -e "$path_destination_parent/error.txt" \
  $path_scripts/4_run_batch_organize_gwas_ldsc_heritability.sh \
  $path_batch_instances \
  $batch_instances_count \
  $path_destination_parent \
  $path_script_gwas_organization \
  $path_scripts \
  $path_promiscuity_scripts
fi
