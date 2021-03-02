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
path_destination_parent=$2 # full path to destination directory
name_prefix=$3 # file name prefix before metabolite identifier or empty string
name_suffix=$4 # file name suffix after metabolite identifier or empty string
file_pattern=$5 # glob pattern by which to recognize relevant files in source directory
path_script_gwas_organization=$6 # full path to script to use for format organization
path_scripts=$7 # full path to scripts for current implementation pipeline
path_promiscuity_scripts=$8 # full path to scripts from promiscuity package

path_batch_instances="${path_destination_parent}/batch_instances.txt"
# Define glob pattern for file paths.
# This definition expands to an array of all files in the path directory that
# match the pattern.
path_pattern="${path_source}/${file_pattern}"

# Initialize directories.
rm -r $path_destination_parent
if [ ! -d $path_destination_parent ]; then
    # Directory does not already exist.
    # Create directory.
    mkdir -p $path_destination_parent
fi

# Report.
echo "----------------------------------------------------------------------"
echo "Report."
echo "----------------------------------------------------------------------"
echo "----------"
echo "source path: " $path_source
#echo "file path pattern: " $path_pattern
echo "destination path: " $path_destination_parent
echo "path to batch instances: " $path_batch_instances
echo "----------"

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
# Iterate on all files and directories in parent directory.
for path_file in $path_source/*; do
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
echo "first batch instance: " ${batch_instances[0]}
echo "last batch instance: " ${batch_instances[batch_instances_count - 1]}

path_file=${batch_instances[0]}
# Determine file name.
file_name="$(basename -- $path_file)"
echo "file: " $file_name

# Determine metabolite identifier.
# Refer to documnetation for test: https://www.freebsd.org/cgi/man.cgi?test
metabolite=${file_name}
if [[ ! -z "$name_prefix" ]]; then
  metabolite=${metabolite/$name_prefix/""}
fi
if [[ ! -z "$name_suffix" ]]; then
  metabolite=${metabolite/$name_suffix/""}
fi
echo "metabolite: " $metabolite





# Submit array batch to Sun Grid Engine.
# Array batch indices cannot start at zero.
# Array batch indices start at one.
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
  $name_prefix \
  $name_suffix \
  $path_script_gwas_organization \
  $path_scripts \
  $path_promiscuity_scripts
fi
